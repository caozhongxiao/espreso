
#include "vtklegacy.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/utilities/utils.h"

#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/output.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/elements/element.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/fetidatastore.h"
#include "../../../mesh/store/surfacestore.h"
#include <fstream>
#include <algorithm>


using namespace espreso;

VTKLegacyDebugInfo::VTKLegacyDebugInfo(const Mesh &mesh, double clusterShrinkRatio, double domainShrinkRatio)
: VTKLegacy(mesh, clusterShrinkRatio, domainShrinkRatio)
{
	_path = Esutils::createDirectory({ Logging::outputRoot(), "VTKLEGACY_DEBUG_OUTPUT" });
}

VTKLegacy::VTKLegacy(const Mesh &mesh, double clusterShrinkRatio, double domainShrinkRatio)
: Visualization(mesh), _clusterShrinkRatio(clusterShrinkRatio), _domainShrinkRatio(domainShrinkRatio)
{

}

void VTKLegacy::mesh(const std::string &name)
{
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	size_t nsize = 0;
	std::vector<eslocal> doffset(_mesh.elements->ndomains);

	for (size_t d = 0; d < _mesh.elements->ndomains; d++) {
		doffset[d] = nsize;
		nsize += _mesh.nodes->dintervals[d].back().DOFOffset + _mesh.nodes->dintervals[d].back().end - _mesh.nodes->dintervals[d].back().begin;
	}

	os << "POINTS " << nsize << " float\n";
	for (size_t d = 0; d < _mesh.elements->ndomains; d++) {
		const auto &coordinates = _mesh.nodes->coordinates->datatarray();
		Point p;
		for (size_t i = 0; i < _mesh.nodes->dintervals[d].size(); ++i) {
			for (eslocal n = _mesh.nodes->dintervals[d][i].begin; n < _mesh.nodes->dintervals[d][i].end; ++n) {
				p = shrink(coordinates[n], _mesh.nodes->center, _mesh.nodes->dcenter[d], _clusterShrinkRatio, _domainShrinkRatio);
				os << p.x << " " << p.y << " " << p.z << "\n";
			}
		}
	}
	os << "\n";

	os << "CELLS " << _mesh.elements->size << " " << _mesh.elements->size + _mesh.elements->nodes->datatarray().size() << "\n";
	auto enodes = _mesh.elements->nodes->cbegin();
	for (eslocal d = 0; d < _mesh.elements->ndomains; d++) {
		for (eslocal e = _mesh.elements->elementsDistribution[d]; e < _mesh.elements->elementsDistribution[d + 1]; ++e, ++enodes) {
			os << enodes->size() << " ";
			for (auto n = enodes->begin(); n != enodes->end(); ++n) {
				auto it = std::lower_bound(_mesh.nodes->dintervals[d].begin(), _mesh.nodes->dintervals[d].end(), *n, [] (const DomainInterval &interval, eslocal node) { return interval.end <= node; });
				os << doffset[d] + it->DOFOffset + *n - it->begin << " ";
			}
			os << "\n";
		}
	}
	os << "\n";

	os << "CELL_TYPES " << _mesh.elements->size << "\n";
	for (auto e = _mesh.elements->epointers->datatarray().begin(); e != _mesh.elements->epointers->datatarray().end(); ++e) {
		switch ((*e)->code) {
		case Element::CODE::SQUARE4:
			os << "9\n";
			break;
		case Element::CODE::SQUARE8:
			os << "23\n";
			break;
		case Element::CODE::TRIANGLE3:
			os << "5\n";
			break;
		case Element::CODE::TRIANGLE6:
			os << "22\n";
			break;
		case Element::CODE::TETRA4:
			os << "10\n";
			break;
		case Element::CODE::TETRA10:
			os << "24\n";
			break;
		case Element::CODE::PYRAMID5:
			os << "14\n";
			break;
		case Element::CODE::PYRAMID13:
			os << "27\n";
			break;
		case Element::CODE::PRISMA6:
			os << "13\n";
			break;
		case Element::CODE::PRISMA15:
			os << "26\n";
			break;
		case Element::CODE::HEXA8:
			os << "12\n";
			break;
		case Element::CODE::HEXA20:
			os << "25\n";
			break;
		default:
			break;
		}
	}
	os << "\n";

	os << "CELL_DATA " << _mesh.elements->size << "\n";
	os << "SCALARS cluster int 1\n";
	os << "LOOKUP_TABLE default\n";
	for (size_t e = 0; e < _mesh.elements->size; e++) {
		os << environment->MPIrank << "\n";
	}
	os << "\n";
	os << "\n";

	os << "SCALARS domains int 1\n";
	os << "LOOKUP_TABLE default\n";
	for (eslocal d = 0; d < _mesh.elements->ndomains; d++) {
		for (eslocal e = _mesh.elements->elementsDistribution[d]; e < _mesh.elements->elementsDistribution[d + 1]; ++e) {
			os << _mesh.elements->firstDomain + d << "\n";
		}
	}
	os << "\n";
}

void VTKLegacy::solution(const std::string &name)
{
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	size_t nsize = 0;
	std::vector<eslocal> doffset(_mesh.elements->ndomains);

	for (size_t d = 0; d < _mesh.elements->ndomains; d++) {
		doffset[d] = nsize;
		nsize += _mesh.nodes->dintervals[d].back().DOFOffset + _mesh.nodes->dintervals[d].back().end - _mesh.nodes->dintervals[d].back().begin;
	}

	os << "POINTS " << nsize << " float\n";
	for (size_t d = 0; d < _mesh.elements->ndomains; d++) {
		const auto &coordinates = _mesh.nodes->coordinates->datatarray();
		Point p;
		for (size_t i = 0; i < _mesh.nodes->dintervals[d].size(); ++i) {
			for (eslocal n = _mesh.nodes->dintervals[d][i].begin; n < _mesh.nodes->dintervals[d][i].end; ++n) {
				p = shrink(coordinates[n], _mesh.nodes->center, _mesh.nodes->dcenter[d], _clusterShrinkRatio, _domainShrinkRatio);
				os << p.x << " " << p.y << " " << p.z << "\n";
			}
		}
	}
	os << "\n";

	os << "CELLS " << _mesh.elements->size << " " << _mesh.elements->size + _mesh.elements->nodes->datatarray().size() << "\n";
	auto enodes = _mesh.elements->nodes->cbegin();
	for (eslocal d = 0; d < _mesh.elements->ndomains; d++) {
		for (eslocal e = _mesh.elements->elementsDistribution[d]; e < _mesh.elements->elementsDistribution[d + 1]; ++e, ++enodes) {
			os << enodes->size() << " ";
			for (auto n = enodes->begin(); n != enodes->end(); ++n) {
				auto it = std::lower_bound(_mesh.nodes->dintervals[d].begin(), _mesh.nodes->dintervals[d].end(), *n, [] (const DomainInterval &interval, eslocal node) { return interval.end <= node; });
				os << doffset[d] + it->DOFOffset + *n - it->begin << " ";
			}
			os << "\n";
		}
	}
	os << "\n";

	os << "CELL_TYPES " << _mesh.elements->size << "\n";
	for (auto e = _mesh.elements->epointers->datatarray().begin(); e != _mesh.elements->epointers->datatarray().end(); ++e) {
		switch ((*e)->code) {
		case Element::CODE::SQUARE4:
			os << "9\n";
			break;
		case Element::CODE::SQUARE8:
			os << "23\n";
			break;
		case Element::CODE::TRIANGLE3:
			os << "5\n";
			break;
		case Element::CODE::TRIANGLE6:
			os << "22\n";
			break;
		case Element::CODE::TETRA4:
			os << "10\n";
			break;
		case Element::CODE::TETRA10:
			os << "24\n";
			break;
		case Element::CODE::PYRAMID5:
			os << "14\n";
			break;
		case Element::CODE::PYRAMID13:
			os << "27\n";
			break;
		case Element::CODE::PRISMA6:
			os << "13\n";
			break;
		case Element::CODE::PRISMA15:
			os << "26\n";
			break;
		case Element::CODE::HEXA8:
			os << "12\n";
			break;
		case Element::CODE::HEXA20:
			os << "25\n";
			break;
		default:
			break;
		}
	}
	os << "\n";

	os << "POINT_DATA " << nsize << "\n";
	os << "SCALARS TEMPERATURE double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (size_t d = 0; d < _mesh.elements->ndomains; d++) {
		for (size_t i = 0; i < _mesh.nodes->dintervals[d].size(); ++i) {
			for (size_t n = 0; n < _mesh.nodes->dintervals[d][i].end - _mesh.nodes->dintervals[d][i].begin; n++) {
				os << (*_mesh.nodes->data.front()->decomposedData)[d][_mesh.nodes->dintervals[d][i].DOFOffset + n] << "\n";
			}
		}
	}
}

void VTKLegacy::nodesIntervals(const std::string &name)
{
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	os << "POINTS " << _mesh.nodes->size << " float\n";
	for (auto n = _mesh.nodes->coordinates->datatarray().begin(); n != _mesh.nodes->coordinates->datatarray().end(); ++n) {
		os << n->x << " " << n->y << " " << n->z << "\n";
	}
	os << "\n";

	os << "CELLS " << _mesh.nodes->size << " " << 2 * _mesh.nodes->size << "\n";
	for (size_t n = 0; n < _mesh.nodes->size; ++n) {
		os << "1 " << n << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << _mesh.nodes->size << "\n";
	for (size_t n = 0; n < _mesh.nodes->size; ++n) {
		os << "1\n";
	}
	os << "\n";

	os << "CELL_DATA " << _mesh.nodes->size << "\n";
	os << "SCALARS interval int 1\n";
	os << "LOOKUP_TABLE default\n";
	int interval = 0;
	for (eslocal i = 0; i < _mesh.nodes->pintervals.size(); i++) {
		for (eslocal n = _mesh.nodes->pintervals[i].begin; n < _mesh.nodes->pintervals[i].end; n++) {
			os << interval << "\n";
		}
		interval++;
	}
	os << "SCALARS sourceProcess int 1\n";
	os << "LOOKUP_TABLE default\n";
	for (eslocal i = 0; i < _mesh.nodes->pintervals.size(); i++) {
		for (eslocal n = _mesh.nodes->pintervals[i].begin; n < _mesh.nodes->pintervals[i].end; n++) {
			os << _mesh.nodes->pintervals[i].sourceProcess << "\n";
		}
	}
	os << "\n";
}

void VTKLegacy::externalIntervals(const std::string &name)
{
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	os << "POINTS " << _mesh.nodes->size << " float\n";
	for (auto n = _mesh.nodes->coordinates->datatarray().begin(); n != _mesh.nodes->coordinates->datatarray().end(); ++n) {
		os << n->x << " " << n->y << " " << n->z << "\n";
	}
	os << "\n";

	eslocal csize = 0;
	for (eslocal i = 0; i < _mesh.nodes->externalIntervals.size(); i++) {
		eslocal ii = _mesh.nodes->externalIntervals[i];
		csize += _mesh.nodes->pintervals[ii].end - _mesh.nodes->pintervals[ii].begin;
	}

	os << "CELLS " << csize << " " << 2 * csize << "\n";
	for (eslocal i = 0; i < _mesh.nodes->externalIntervals.size(); i++) {
		eslocal ii = _mesh.nodes->externalIntervals[i];
		for (eslocal n = _mesh.nodes->pintervals[ii].begin; n < _mesh.nodes->pintervals[ii].end; n++) {
			os << "1 " << n << "\n";
		}
	}
	os << "\n";

	os << "CELL_TYPES " << csize << "\n";
	for (size_t n = 0; n < csize; ++n) {
		os << "1\n";
	}
	os << "\n";

	os << "CELL_DATA " << csize << "\n";
	os << "SCALARS interval int 1\n";
	os << "LOOKUP_TABLE default\n";
	int interval = 0;
	for (eslocal i = 0; i < _mesh.nodes->externalIntervals.size(); i++) {
		eslocal ii = _mesh.nodes->externalIntervals[i];
		for (eslocal n = _mesh.nodes->pintervals[ii].begin; n < _mesh.nodes->pintervals[ii].end; n++) {
			os << interval << "\n";
		}
		interval++;
	}
	os << "\n";
}

void VTKLegacy::sharedInterface(const std::string &name)
{
	if (_mesh.FETIData == NULL || _mesh.FETIData->interfaceNodes == NULL) {
		return;
	}
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	size_t points = _mesh.FETIData->interfaceNodes->datatarray().size();

	os << "POINTS " << points << " float\n";
	for (size_t i = 0; i < points; i++) {
		const Point &n = _mesh.nodes->coordinates->datatarray()[_mesh.FETIData->interfaceNodes->datatarray()[i]];
		os << n.x << " " << n.y << " " << n.z << "\n";
	}
	os << "\n";

	os << "CELLS " << points << " " << 2 * points << "\n";
	for (size_t n = 0; n < points; ++n) {
		os << "1 " << n << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << points << "\n";
	for (size_t n = 0; n < points; ++n) {
		os << "1\n";
	}
	os << "\n";

	os << "CELL_DATA " << points << "\n";
	os << "SCALARS distribution int 1\n";
	os << "LOOKUP_TABLE default\n";
	for (eslocal i = 0; i < _mesh.FETIData->inodesDistribution.size() - 1; i++) {
		for (eslocal n = _mesh.FETIData->inodesDistribution[i]; n < _mesh.FETIData->inodesDistribution[i + 1]; n++) {
			os << i << "\n";
		}
	}
	os << "\n";
}

void VTKLegacy::domainSurface(const std::string &name)
{
	if (_mesh.domainsSurface == NULL) {
		return;
	}

	auto surface = [&] (const std::string &suffix, serializededata<eslocal, eslocal>* elements, std::vector<eslocal> &distribution) {
		std::ofstream os(name + "." + suffix + std::to_string(environment->MPIrank) + ".vtk");

		os << "# vtk DataFile Version 2.0\n";
		os << "EXAMPLE\n";
		os << "ASCII\n";
		os << "DATASET UNSTRUCTURED_GRID\n\n";

		os << "POINTS " << _mesh.domainsSurface->coordinates->datatarray().size() << " float\n";
		for (size_t d = 0; d < _mesh.elements->ndomains; ++d) {
			for (eslocal n = _mesh.domainsSurface->cdistribution[d]; n < _mesh.domainsSurface->cdistribution[d + 1]; ++n) {
				Point p = shrink(_mesh.domainsSurface->coordinates->datatarray()[n], _mesh.nodes->center, _mesh.nodes->dcenter[d], _clusterShrinkRatio, _domainShrinkRatio);
				os << p.x << " " << p.y << " " << p.z << "\n";
			}
		}
		os << "\n";

		os << "CELLS " << elements->structures() << " " << elements->structures() + elements->datatarray().size() << "\n";
		for (size_t d = 0; d < _mesh.elements->ndomains; ++d) {
			for (auto element = elements->cbegin() + distribution[d]; element != elements->cbegin() + distribution[d + 1]; ++element) {
				os << element->size() << " ";
				for (auto n = element->begin(); n != element->end(); ++n) {
					os << _mesh.domainsSurface->cdistribution[d] + *n << " ";
				}
				os << "\n";
			}
		}
		os << "\n";

		os << "CELL_TYPES " << distribution.back() << "\n";
		for (size_t n = 0; n < distribution.back(); ++n) {
			os << "5\n";
		}
		os << "\n";
	};

	surface("elements",  _mesh.domainsSurface->elements, _mesh.domainsSurface->edistribution);
	if (_mesh.domainsSurface->triangles != NULL) {
		surface("triangles",  _mesh.domainsSurface->triangles, _mesh.domainsSurface->tdistribution);
	}
}

void VTKLegacy::points(const std::string &name, const std::vector<eslocal> &points, const std::vector<eslocal> &distribution)
{
	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	os << "POINTS " << points.size() << " float\n";
	for (size_t d = 0; d < distribution.size() - 1; d++) {
		for (size_t i = distribution[d]; i < distribution[d + 1]; i++) {
			Point n = shrink(_mesh.nodes->coordinates->datatarray()[points[i]], _mesh.nodes->center, _mesh.nodes->dcenter[d], _clusterShrinkRatio, _domainShrinkRatio);
			os << n.x << " " << n.y << " " << n.z << "\n";
		}
	}
	os << "\n";

	os << "CELLS " << points.size() << " " << 2 * points.size() << "\n";
	for (size_t n = 0; n < points.size(); ++n) {
		os << "1 " << n << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << points.size() << "\n";
	for (size_t n = 0; n < points.size(); ++n) {
		os << "1\n";
	}
	os << "\n";
}

void VTKLegacy::corners(const std::string &name)
{
	if (_mesh.FETIData == NULL) {
		return;
	}

	std::ofstream os(name + std::to_string(environment->MPIrank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	eslocal points = 0;
	auto cdomains = _mesh.FETIData->cornerDomains->begin();
	for (size_t i = 0; i < _mesh.FETIData->corners.size(); ++i, ++cdomains) {
		points += cdomains->size();
	}
	os << "POINTS " << points << " float\n";
	cdomains = _mesh.FETIData->cornerDomains->begin();
	for (size_t i = 0; i < _mesh.FETIData->corners.size(); ++i, ++cdomains) {
		for (auto d = cdomains->begin(); d != cdomains->end(); ++d) {
			Point n = shrink(_mesh.nodes->coordinates->datatarray()[_mesh.FETIData->corners[i]], _mesh.nodes->center, _mesh.nodes->dcenter[*d], _clusterShrinkRatio, _domainShrinkRatio);
			os << n.x << " " << n.y << " " << n.z << "\n";
		}
	}
	os << "\n";

	os << "CELLS " << points << " " << 2 * points << "\n";
	for (size_t n = 0; n < points; ++n) {
		os << "1 " << n << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << points << "\n";
	for (size_t n = 0; n < points; ++n) {
		os << "1\n";
	}
	os << "\n";
}

void VTKLegacy::sFixPoints(const std::string &name)
{
	if (_mesh.FETIData == NULL) {
		return;
	}
	points(name, _mesh.FETIData->surfaceFixPoints, _mesh.FETIData->sFixPointsDistribution);
}

void VTKLegacy::iFixPoints(const std::string &name)
{
	if (_mesh.FETIData == NULL) {
		return;
	}
	points(name, _mesh.FETIData->innerFixPoints, _mesh.FETIData->iFixPointsDistribution);
}

