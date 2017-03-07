
#include "resultstore.h"

#include <numeric>

#include "../configuration/output.h"
#include "../configuration/environment.h"

#include "../basis/point/point.h"
#include "../mesh/elements/element.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"
#include "../mesh/settings/property.h"
#include "../mesh/structures/elementtypes.h"

#include "../assembler/step.h"
#include "../assembler/solution.h"

#include "../basis/utilities/utils.h"

using namespace espreso::output;

DataArrays::~DataArrays()
{
	for (auto it = elementDataDouble.begin(); it != elementDataDouble.end(); ++it) {
		delete it->second;
	}
	for (auto it = elementDataInteger.begin(); it != elementDataInteger.end(); ++it) {
		delete it->second;
	}
	for (auto it = pointDataDouble.begin(); it != pointDataDouble.end(); ++it) {
		delete it->second;
	}
	for (auto it = pointDataInteger.begin(); it != pointDataInteger.end(); ++it) {
		delete it->second;
	}
}

ResultStore::ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: Store(output), _mesh(mesh), _path(path)
{

}

ResultStore::~ResultStore()
{

}

static espreso::Point computeClusterCenter(const espreso::Mesh *mesh)
{
	espreso::Point clusterCenter;
	for (size_t i = 0; i < mesh->coordinates().clusterSize(); i++) {
		clusterCenter += mesh->coordinates()[i];
	}
	clusterCenter /= mesh->coordinates().clusterSize();
	return clusterCenter;
}

static std::vector<espreso::Point> computeDomainsCenters(const espreso::Mesh *mesh)
{
	std::vector<espreso::Point> domainsCenters(mesh->parts());
	for (size_t p = 0; p < mesh->coordinates().parts(); p++) {
		for (size_t i = 0; i < mesh->coordinates().localSize(p); i++) {
			domainsCenters[p] += mesh->coordinates().get(i, p);
		}
		domainsCenters[p] /= mesh->coordinates().localSize(p);
	}
	return domainsCenters;
}

void ResultStore::coordinatePreprocessing(const std::vector<std::vector<eslocal> > &indices, std::vector<double> &coordinates, std::vector<size_t> &offsets)
{
	// compute coordinates
	std::vector<std::vector<double> > dCoordinates(_mesh->parts());
	#pragma omp parallel for
	for (size_t p = 0; p < indices.size(); p++) {
		Point point;
		dCoordinates[p].reserve(indices[p].size());
		for (size_t i = 0; i < indices[p].size(); i++) {
			point = _mesh->coordinates()[indices[p][i]];
			point = _clusterCenter + (point - _clusterCenter) * _configuration.cluster_shrink_ratio;
			point = _domainsCenters[p] + (point - _domainsCenters[p]) * _configuration.domain_shrink_ratio;
			dCoordinates[p].insert(dCoordinates[p].end(), { point.x, point.y, point.z });
		}
	}

	coordinates.clear();
	offsets = { 0 };
	for (size_t p = 0; p < _mesh->parts(); p++) {
		coordinates.insert(coordinates.end(), dCoordinates[p].begin(), dCoordinates[p].end());
		offsets.push_back(coordinates.size() / 3);
	}
}

void ResultStore::elementsPreprocessing(const std::vector<Element*> &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	if (!region.size()) {
		// TODO: collected output
		return;
	}

	std::vector<std::vector<eslocal> > rCoordinates(_mesh->parts());
	for (size_t e = 0; e < region.size(); e++) {
		for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d) {
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				rCoordinates[*d].push_back(region[e]->node(n));
			}
		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		std::sort(rCoordinates[p].begin(), rCoordinates[p].end());
		Esutils::removeDuplicity(rCoordinates[p]);
	}

	std::vector<size_t> offsets;
	coordinatePreprocessing(rCoordinates, coordinates, offsets);

	for (size_t e = 0, offset = 0; e < region.size(); e++) {
		elementsTypes.insert(elementsTypes.end(), region[e]->domains().size(), region[e]->vtkCode());
		for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d, offset += region[e]->nodes()) {
			elementsNodes.push_back(offset + region[e]->nodes());
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				eslocal oIndex = std::lower_bound(rCoordinates[*d].begin(), rCoordinates[*d].end(), region[e]->node(n)) - rCoordinates[*d].begin();
				elements.push_back(oIndex + offsets[*d]);
			}
		}
	}
}

void ResultStore::regionData(size_t step, const espreso::Region &region, DataArrays &data)
{
	if (!region.elements().size()) {
		return;
	}
	for (auto it = region.settings[step].begin(); it != region.settings[step].end(); ++it) {
		std::vector<double> *values = new std::vector<double>();
		values->reserve(region.elements().size());
		for (size_t e = 0; e < region.elements().size(); e++) {
			values->push_back(0);
			for (size_t n = 0; n < region.elements()[e]->nodes(); n++) {
				values->back() += region.elements()[e]->sumProperty(it->first, n, step, 0);
			}
			values->back() /= region.elements()[e]->nodes();
			values->insert(values->end(), region.elements()[e]->domains().size() - 1, values->back());
		}
		std::stringstream ss;
		ss << it->first;
		data.elementDataDouble[ss.str()] = values;
	}
}

void ResultStore::preprocessing()
{
	if (_configuration.collected) {
		// Implement collected printing
	} else {
		_clusterCenter = computeClusterCenter(_mesh);
		_domainsCenters = computeDomainsCenters(_mesh);
		std::vector<size_t> domainOffset;

		coordinatePreprocessing(_mesh->coordinates().localToCluster(), _coordinates, domainOffset);

		// compute elements
		std::vector<std::vector<eslocal> > dElementsTypes(_mesh->parts());
		std::vector<std::vector<eslocal> > dElementsNodes(_mesh->parts());
		std::vector<std::vector<eslocal> > dElements(_mesh->parts());
		std::vector<size_t> offsets(_mesh->parts());
		#pragma omp parallel for
		for (size_t p = 0; p < _mesh->parts(); p++) {
			size_t offset = 0;
			for (size_t e = _mesh->getPartition()[p]; e < _mesh->getPartition()[p + 1]; e++) {
				dElementsTypes[p].push_back(_mesh->elements()[e]->vtkCode());
				dElementsNodes[p].push_back(_mesh->elements()[e]->nodes());
				offset += _mesh->elements()[e]->nodes();
				for (size_t n = 0; n < _mesh->elements()[e]->nodes(); n++) {
					dElements[p].push_back(_mesh->coordinates().localIndex(_mesh->elements()[e]->node(n), p) + domainOffset[p]);
				}
			}
			offsets[p] = offset;
		}

		Esutils::sizesToOffsets(offsets);
		#pragma omp parallel for
		for (size_t p = 0; p < _mesh->parts(); p++) {
			for (size_t e = _mesh->getPartition()[p], i = 0; e < _mesh->getPartition()[p + 1]; e++, i++) {
				dElementsNodes[p][i] = offsets[p] += dElementsNodes[p][i];
			}
		}

		for (size_t p = 0; p < _mesh->parts(); p++) {
			_elementsTypes.insert(_elementsTypes.end(), dElementsTypes[p].begin(), dElementsTypes[p].end());
			_elementsNodes.insert(_elementsNodes.end(), dElementsNodes[p].begin(), dElementsNodes[p].end());
			_elements.insert(_elements.end(), dElements[p].begin(), dElements[p].end());
		}
	}
}

static void fillMeshSettings(DataArrays &data, const std::vector<double> &coordinates, const std::vector<eslocal> &partPtrs)
{
	std::vector<eslocal> *pointID = new std::vector<eslocal>(coordinates.size() / 3);
	std::vector<eslocal> *elementID = new std::vector<eslocal>(partPtrs.back());
	std::vector<eslocal> *decomposition = new std::vector<eslocal>();

	std::iota(pointID->begin(), pointID->end(), 0);
	std::iota(elementID->begin(), elementID->end(), 0);
	for (size_t p = 0; p < partPtrs.size() - 1; p++) {
		decomposition->insert(decomposition->end(), partPtrs[p + 1] - partPtrs[p], p);
	}

	data.pointDataInteger["pointID"] = pointID;
	data.elementDataInteger["elementID"] = elementID;
	data.elementDataInteger["decomposition"] = decomposition;
}

void ResultStore::storeSettings(const Step &step)
{
	storeSettings(std::vector<size_t>{ step.step });
}

void ResultStore::storeSettings(size_t steps)
{
	std::vector<size_t> _steps(steps);
	std::iota(_steps.begin(), _steps.end(), 0);
	storeSettings(_steps);
}

void ResultStore::storeSettings(const std::vector<size_t> &steps)
{
	std::vector<std::string> prefixes;
	std::vector<std::string> roots;

	for (size_t step = 0; step < steps.size(); step++) {
		if (!environment->MPIrank) {
			roots.push_back(Esutils::createDirectory({ "results", "step" + std::to_string(steps[step]) }));
		}
		prefixes.push_back(Esutils::createDirectory({ "results", "step" + std::to_string(steps[step]), std::to_string(environment->MPIrank) }));
	}

	DataArrays data;
	fillMeshSettings(data, _coordinates, _mesh->getPartition());
	for (size_t step = 0; step < steps.size(); step++) {
		store(prefixes[step] + "mesh", _coordinates, _elementsTypes, _elementsNodes, _elements, data);
		if (!environment->MPIrank) {
			linkClusters(roots[step], "mesh", data);
		}
	}

	for (size_t r = 2; r < _mesh->regions().size(); r++) {

		std::vector<double> coordinates;
		std::vector<eslocal> elementTypes, elementNodes, elements;
		elementsPreprocessing(_mesh->regions()[r]->elements(), coordinates, elementTypes, elementNodes, elements);

		for (size_t step = 0; step < steps.size(); step++) {
			DataArrays rData;
			regionData(steps[step], *_mesh->regions()[r], rData);

			store(prefixes[step] + _mesh->regions()[r]->name, coordinates, elementTypes, elementNodes, elements, rData);
			if (!environment->MPIrank) {
				linkClusters(roots[step], _mesh->regions()[r]->name, rData);
			}
		}
	}
}

void ResultStore::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	Step step;
	std::vector<Solution*> solution;
	std::vector<Property> props;
	if (dimension == 1) {
		props.push_back(Property::TEMPERATURE);
	}
	if (dimension == 3) {
		props.push_back(Property::DISPLACEMENT_X);
		props.push_back(Property::DISPLACEMENT_Y);
		props.push_back(Property::DISPLACEMENT_Z);
	}

	solution.push_back(new Solution(name, eType, props, values));
	storeSolution(step, solution);
	delete solution.back();
	if (!environment->MPIrank) {
		linkSteps("solution", _steps);
	}
	_steps.clear();
}

void ResultStore::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	std::string prefix;
	std::string root;

	if (!environment->MPIrank) {
		root = Esutils::createDirectory({ "results", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep) });
	}
	prefix = Esutils::createDirectory({ "results", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), std::to_string(environment->MPIrank) });

	store(prefix + "solution", _coordinates, _elementsTypes, _elementsNodes, _elements, solution);
	if (!environment->MPIrank) {
		linkClusters(root, "solution", solution, _coordinates.size() / 3, _elementsTypes.size());
	}
	_steps.push_back(std::make_pair(root, step));
}

void ResultStore::finalize()
{
	if (!environment->MPIrank && _steps.size()) {
		linkSteps("solution", _steps);
	}
}


void ResultStore::storeFETIData(const Step &step, const Instance &instance)
{
	storeFixPoints(step);
	storeCorners(step);
	storeDirichlet(step, instance);
	storeLambdas(step, instance);
}

void ResultStore::storeFixPoints(const Step &step)
{

}

void ResultStore::storeCorners(const Step &step)
{

}

void ResultStore::storeDirichlet(const Step &step, const Instance &instance)
{

}

void ResultStore::storeLambdas(const Step &step, const Instance &instance)
{

}






