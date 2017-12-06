
#include "ensight.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/logging/logging.h"
#include "../../../basis/utilities/communication.h"
#include "../../../basis/utilities/utils.h"
#include "../../../basis/utilities/parser.h"

#include "../../../config/ecf/environment.h"

#include "../../../mesh/elements/element.h"
#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"
#include "../../../mesh/store/boundaryregionstore.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <functional>

using namespace espreso;

EnSight::EnSight(const std::string &name, const Mesh &mesh)
: _path(Logging::outputRoot() + "/"), _name(name), _mesh(mesh), _casefile(NULL)
{
	if (environment->MPIrank == 0) {
		_casefile = new std::ofstream(_path + _name + ".case");
		(*_casefile) << "#\n";
		(*_casefile) << "# ESPRESO solution\n";
		(*_casefile) << "#\n";

		(*_casefile) << "\nFORMAT\n";
		(*_casefile) << "type: \tensight gold\n\n";

		(*_casefile) << "GEOMETRY\n\n";
		(*_casefile) << "model:\t" << _name << ".geo\n\n";

		(*_casefile) << "VARIABLE\n\n";

		(*_casefile).flush();
	}
}

EnSight::~EnSight()
{
	if (_casefile != NULL) {
		delete _casefile;
	}
}

std::string EnSight::codetotype(int code)
{
	switch (static_cast<Element::CODE>(code)) {

	case Element::CODE::POINT1: return "point";

	case Element::CODE::LINE2: return "bar2";

	case Element::CODE::TRIANGLE3: return "tria3";
	case Element::CODE::SQUARE4: return "quad4";

	case Element::CODE::TETRA4: return "tetra4";
	case Element::CODE::PYRAMID5: return "pyramid5";
	case Element::CODE::PRISMA6: return "penta6";
	case Element::CODE::HEXA8: return "hexa8";

	case Element::CODE::LINE3: return "bar3";

	case Element::CODE::TRIANGLE6: return "tria6";
	case Element::CODE::SQUARE8: return "quad8";

	case Element::CODE::TETRA10: return "tetra10";
	case Element::CODE::PYRAMID13: return "pyramid13";
	case Element::CODE::PRISMA15: return "penta15";
	case Element::CODE::HEXA20: return "hexa20";

	default:
		ESINFO(ERROR) << "ESPRESO internal error: unknown element code.";
		return "";
	}
}

void EnSight::storeGeometry()
{
	std::string name = _path + _name + ".geo";
	int part = 1;

	std::stringstream os;
	if (environment->MPIrank == 0) {
		os << "EnSight Gold geometry format\n";
		os << "----------------------------\n";

		os << "node id off\n";
		os << "element id off\n";
	}

	auto storePartHeader = [&] (const std::string &name, eslocal nodes) {
		if (environment->MPIrank == 0) {
			os << "part\n";
			os << std::setw(10) << part++ << "\n";
			os << name << "\n";

			os << "coordinates\n";
			os << std::setw(10) << nodes << "\n";
		}
	};

	auto storeRegionNodes = [&] (const std::vector<ProcessInterval> &intervals, const serializededata<eslocal, eslocal> *nodes, std::function<double(const Point *p)> getCoordinate) {
		os << std::showpos << std::scientific << std::setprecision(5);
		for (size_t i = 0; i < intervals.size(); i++) {
			if (intervals[i].sourceProcess == environment->MPIrank) {
				for (auto n = nodes->datatarray().cbegin() + intervals[i].begin; n != nodes->datatarray().cbegin() + intervals[i].end; ++n) {
					os << getCoordinate(_mesh.nodes->coordinates->datatarray().cbegin() + *n) << "\n";
				}
			}
		}
		os << std::noshowpos;
		pushInterval(os.str().size());
	};

	auto storeRegionElements = [&] (
			const std::vector<eslocal> &ecounters,
			const tarray<eslocal> &elements, const std::vector<ElementsInterval> &eintervals,
			const tarray<eslocal> &nodes, const std::vector<ProcessInterval> &nintervals) {

		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (ecounters[etype]) {
				if (environment->MPIrank == 0) {
					os << codetotype(etype) << "\n";
					os << std::setw(10) << ecounters[etype] << "\n";
				}

				for (size_t i = 0; i < eintervals.size(); i++) {
					if (eintervals[i].code == etype) {
						auto enodes = _mesh.elements->nodes->cbegin() + eintervals[i].begin;
						eslocal prev = eintervals[i].begin;
						for (eslocal e = eintervals[i].begin; e < eintervals[i].end; prev = elements[e++]) {
							enodes += elements[e] - prev;
							for (auto n = enodes->begin(); n != enodes->end(); ++n) {
								auto iit = std::lower_bound(_mesh.nodes->pintervals.begin(), _mesh.nodes->pintervals.end(), *n, [] (const ProcessInterval &interval, eslocal node) { return interval.end <= node; });
								size_t iindex = iit - _mesh.nodes->pintervals.begin();
								eslocal offset = std::lower_bound(nodes.begin() + nintervals[iindex].begin, nodes.begin() + nintervals[iindex].end, *n) - nodes.begin();
								os << std::setw(10) << nintervals[iindex].globalOffset + offset - nintervals[iindex].begin + 1;
							}
							os << "\n";
						}
					}
				}
				pushInterval(os.str().size());
			}
		}
	};

	auto storeERegion = [&] (const ElementsRegionStore *region) {
		storeRegionNodes(region->nintervals, region->nodes, [] (const Point *p) { return p->x; });
		storeRegionNodes(region->nintervals, region->nodes, [] (const Point *p) { return p->y; });
		storeRegionNodes(region->nintervals, region->nodes, [] (const Point *p) { return p->z; });

		if (StringCompare::caseInsensitiveEq(region->name, "ALL_ELEMENTS")) {
			storeRegionElements(region->ecounters, region->uniqueElements->datatarray(), region->ueintervals, region->nodes->datatarray(), region->nintervals);
		} else {
			storeRegionElements(region->ecounters, region->elements->datatarray(), region->eintervals, region->nodes->datatarray(), region->nintervals);
		}
	};

	auto storeBRegion = [&] () {

	};

	for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
		storePartHeader(_mesh.elementsRegions[r]->name, _mesh.elementsRegions[r]->uniqueTotalSize);
		storeERegion(_mesh.elementsRegions[r]);
	}

	storeIntervals(name, os.str(), commitIntervals());
}

void EnSight::storeFETIData()
{
	auto iterateElements = [&] (std::stringstream &os, const std::vector<ElementsInterval> &intervals, const std::vector<eslocal> &ecounters, std::function<double(eslocal domain)> fnc) {
		os << std::scientific << std::setprecision(5);
		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (ecounters[etype]) {
				if (environment->MPIrank == 0) {
					os << codetotype(etype) << "\n";
				}

				for (size_t i = 0; i < intervals.size(); i++) {
					if (intervals[i].code == etype) {
						auto domain = std::lower_bound(_mesh.elements->elementsDistribution.begin(), _mesh.elements->elementsDistribution.end(), intervals[i].begin) - _mesh.elements->elementsDistribution.begin();
						for (eslocal e = intervals[i].begin; e < intervals[i].end; ++e) {
							os << " " << fnc(domain) << "\n";
						}
					}
				}

				pushInterval(os.str().size());
			}
		}
	};

	eslocal part = 1;

	auto storePartHeader = [&] (std::stringstream &os) {
		if (environment->MPIrank == 0) {
			os << "part\n";
			os << std::setw(10) << part++ << "\n";
		}
	};


	{ // DOMAINS
		std::string filename = _name + ".DOMAINS";
		std::string name = _path + filename;

		std::stringstream os;
		if (environment->MPIrank == 0) {
			(*_casefile) << "scalar per element:\tDOMAINS\t\t" << filename << "\n";
			_casefile->flush();
			os << "DOMAINS\n";
		}

		part = 1;
		clearIntervals();
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader(os);
			if (StringCompare::caseInsensitiveEq(_mesh.elementsRegions[r]->name, "ALL_ELEMENTS")) {
				iterateElements(os, _mesh.elementsRegions[r]->ueintervals, _mesh.elementsRegions[r]->ecounters, [&] (eslocal domain)->double { return domain + _mesh.elements->firstDomain; });
			} else {
				iterateElements(os, _mesh.elementsRegions[r]->eintervals, _mesh.elementsRegions[r]->ecounters, [&] (eslocal domain)->double { return domain + _mesh.elements->firstDomain; });
			}
		}

		storeIntervals(name, os.str(), commitIntervals());
	}

	{ // CLUSTERS
		std::string filename = _name + ".CLUSTERS";
		std::string name = _path + filename;

		std::stringstream os;
		if (environment->MPIrank == 0) {
			(*_casefile) << "scalar per element:\tCLUSTERS\t" << filename << "\n";
			_casefile->flush();
			os << "CLUSTERS\n";
		}

		part = 1;
		clearIntervals();
		eslocal cluster = _mesh.elements->gatherClustersDistribution()[environment->MPIrank];
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader(os);
			if (StringCompare::caseInsensitiveEq(_mesh.elementsRegions[r]->name, "ALL_ELEMENTS")) {
				iterateElements(os, _mesh.elementsRegions[r]->ueintervals, _mesh.elementsRegions[r]->ecounters, [&] (eslocal domain)->double { return _mesh.elements->clusters[domain] + cluster; });
			} else {
				iterateElements(os, _mesh.elementsRegions[r]->eintervals, _mesh.elementsRegions[r]->ecounters, [&] (eslocal domain)->double { return _mesh.elements->clusters[domain] + cluster; });
			}
		}

		storeIntervals(name, os.str(), commitIntervals());
	}
}

void EnSight::storeVariables()
{
	if (_mesh.nodes->data.size() != 1) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement store variables.";
	}

	std::string filename = _name + "." + _mesh.nodes->data.front()->names.front().substr(0, 4);
	std::string name = _path + filename;

	if (environment->MPIrank == 0) {
		(*_casefile) << "scalar per node:\t" << _mesh.nodes->data.front()->names.front() << "\t" << filename << "\n";
		_casefile->flush();
	}

	std::stringstream os;
	if (environment->MPIrank == 0) {
		os << _mesh.nodes->data.front()->names.front() << "\n";
	}

	eslocal part = 1;

	auto storePartHeader = [&] () {
		if (environment->MPIrank == 0) {
			os << "part\n";
			os << std::setw(10) << part++ << "\n";
			os << "coordinates\n";
		}
	};

	auto iterateNodes = [&] (const std::vector<ProcessInterval> &intervals, const tarray<eslocal> &nodes) {
		os << std::showpos << std::scientific << std::setprecision(5);
		for (size_t i = 0; i < intervals.size(); i++) {
			if (intervals[i].sourceProcess == environment->MPIrank) {
				eslocal offset = _mesh.nodes->pintervals[i].globalOffset - _mesh.nodes->uniqueOffset;
				for (eslocal n = intervals[i].begin; n < intervals[i].end; ++n) {
					eslocal index = offset + nodes[n] - _mesh.nodes->pintervals[i].begin;
					os << (*_mesh.nodes->data.front()->gathredData)[index] << "\n";
				}
			}
		}
		os << std::noshowpos;
		pushInterval(os.str().size());
	};

	clearIntervals();
	for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
		storePartHeader();
		iterateNodes(_mesh.elementsRegions[r]->nintervals, _mesh.elementsRegions[r]->nodes->datatarray());
	}

	storeIntervals(name, os.str(), commitIntervals());
}


