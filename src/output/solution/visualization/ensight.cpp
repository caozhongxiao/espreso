
#include "ensight.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/logging/logging.h"
#include "../../../basis/utilities/communication.h"
#include "../../../basis/utilities/utils.h"

#include "../../../config/ecf/environment.h"

#include "../../../mesh/elements/element.h"
#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
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
	std::vector<eslocal> ndistribution = _mesh.nodes->gatherUniqueNodeDistribution();
	std::vector<eslocal> edistribution = _mesh.elements->gatherElementsProcDistribution();

	int part = 1;

	std::stringstream os;
	if (environment->MPIrank == 0) {
		os << "EnSight Gold geometry format\n";
		os << "----------------------------\n";

		os << "node id off\n";
		os << "element id off\n";

		os << "part\n";
		os << std::setw(10) << part++ << "\n";
		os << "MESH\n";

		os << "coordinates\n";
		os << std::setw(10) << ndistribution.back() << "\n";
	}

	auto storeNodes = [&] (std::function<double(const Point *p)> getCoordinate) {
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i) {
			if (_mesh.nodes->pintervals[i].sourceProcess == environment->MPIrank) {
				auto begin = _mesh.nodes->coordinates->datatarray().begin() + _mesh.nodes->pintervals[i].begin;
				auto end = _mesh.nodes->coordinates->datatarray().begin() + _mesh.nodes->pintervals[i].end;
				for (auto n = begin; n != end; ++n) {
					os << std::scientific << std::setprecision(5) << getCoordinate(n) << "\n";
				}
			}
		}
		pushInterval(os.str().size());
	};

	auto storeRegionNodes = [&] (const BoundaryRegionStore *store, std::function<double(const Point *p)> getCoordinate) {
		for (size_t i = 0; i < store->nodesIntervals.size(); i++) {
			if (store->nodesIntervals[i].sourceProcess == environment->MPIrank) {
				for (auto n = store->nodes->datatarray().cbegin() + store->nodesIntervals[i].begin; n != store->nodes->datatarray().cbegin() + store->nodesIntervals[i].end; ++n) {
					os << std::scientific << std::setprecision(5) << getCoordinate(_mesh.nodes->coordinates->datatarray().cbegin() + *n) << "\n";
				}
			}
		}
		pushInterval(os.str().size());
	};

	os << std::showpos;
	storeNodes([] (const Point *p) { return p->x; });
	storeNodes([] (const Point *p) { return p->y; });
	storeNodes([] (const Point *p) { return p->z; });
	os << std::noshowpos;

	for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
		if (_mesh.elements->ecounters[etype]) {
			if (environment->MPIrank == 0) {
				os << codetotype(etype) << "\n";
				os << std::setw(10) << _mesh.elements->ecounters[etype] << "\n";
			}

			for (size_t i = 0; i < _mesh.elements->eintervals.size(); i++) {
				if (_mesh.elements->eintervals[i].code == etype) {
					auto nodes = _mesh.elements->nodes->cbegin() + _mesh.elements->eintervals[i].begin;
					for (eslocal e = _mesh.elements->eintervals[i].begin; e < _mesh.elements->eintervals[i].end; ++e, ++nodes) {
						for (auto n = nodes->begin(); n != nodes->end(); ++n) {
							auto it = std::lower_bound(_mesh.nodes->pintervals.begin(), _mesh.nodes->pintervals.end(), *n, [] (const ProcessInterval &interval, eslocal node) { return interval.end <= node; });
							os << std::setw(10) << it->globalOffset + *n - it->begin + 1;
						}
						os << "\n";
					}
				}
			}

			pushInterval(os.str().size());
		}
	}

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		if (_mesh.boundaryRegions[r]->nodes) {
			if (environment->MPIrank == 0) {
				os << "part\n";
				os << std::setw(10) << part++ << "\n";
				os << _mesh.boundaryRegions[r]->name << "\n";

				os << "coordinates\n";
				os << std::setw(10) << _mesh.boundaryRegions[r]->uniqueTotalSize << "\n";
			}

			os << std::showpos;
			storeRegionNodes(_mesh.boundaryRegions[r], [] (const Point *p) { return p->x; });
			storeRegionNodes(_mesh.boundaryRegions[r], [] (const Point *p) { return p->y; });
			storeRegionNodes(_mesh.boundaryRegions[r], [] (const Point *p) { return p->z; });
			os << std::noshowpos;

			if (environment->MPIrank == 0) {
				os << codetotype(static_cast<int>(Element::CODE::POINT1)) << "\n";
				os << std::setw(10) << _mesh.boundaryRegions[r]->uniqueTotalSize << "\n";
			}

			eslocal offset = 1;
			for (size_t i = 0; i < _mesh.boundaryRegions[r]->nodesIntervals.size(); i++) {
				if (_mesh.boundaryRegions[r]->nodesIntervals[i].sourceProcess == environment->MPIrank) {
					auto begin = _mesh.boundaryRegions[r]->nodes->datatarray().cbegin() + _mesh.boundaryRegions[r]->nodesIntervals[i].begin;
					auto end = _mesh.boundaryRegions[r]->nodes->datatarray().cbegin() + _mesh.boundaryRegions[r]->nodesIntervals[i].end;
					for (auto n = begin; n != end; ++n) {
						os << std::setw(10) << _mesh.boundaryRegions[r]->uniqueOffset + offset++ << "\n";
					}
				}
			}

			pushInterval(os.str().size());
		}
	}

	storeIntervals(name, os.str(), commitIntervals());
}

void EnSight::storeFETIData()
{
	auto iterateElements = [&] (std::stringstream &os, std::function<double(eslocal domain)> fnc) {
		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (_mesh.elements->ecounters[etype]) {
				if (environment->MPIrank == 0) {
					os << codetotype(etype) << "\n";
				}

				for (size_t i = 0; i < _mesh.elements->eintervals.size(); i++) {
					if (_mesh.elements->eintervals[i].code == etype) {
						auto dit = std::lower_bound(_mesh.elements->elementsDistribution.begin(), _mesh.elements->elementsDistribution.end(), _mesh.elements->eintervals[i].begin);
						eslocal d = dit - _mesh.elements->elementsDistribution.begin();
						eslocal begin = *dit;
						if (begin < _mesh.elements->eintervals[i].begin) {
							begin = _mesh.elements->eintervals[i].begin;
							--d;
						}
						if (*(dit + 1) < _mesh.elements->eintervals[i].end) {
							while (begin < _mesh.elements->eintervals[i].end) {
								for (eslocal e = begin; e < *(dit + 1); ++e) {
									os << " " << std::scientific << std::setprecision(5) << fnc(d) << "\n";
								}
								++dit;
								begin = *dit;
								++d;
							}
						} else {
							for (eslocal e = begin; e < _mesh.elements->eintervals[i].end; ++e) {
								os << " " << std::scientific << std::setprecision(5) << fnc(d) << "\n";
							}
						}
					}
				}

				pushInterval(os.str().size());
			}
		}
	};


	{ // DOMAINS
		std::string filename = _name + ".DOMAINS";
		std::string name = _path + filename;

		if (environment->MPIrank == 0) {
			(*_casefile) << "scalar per element:\tDOMAINS\t\t" << filename << "\n";
			_casefile->flush();
		}

		std::stringstream os;

		if (environment->MPIrank == 0) {
			os << "DOMAINS\n";
			os << "part\n";
			os << std::setw(10) << 1 << "\n";
		}

		iterateElements(os, [&] (eslocal domain)->double { return domain + _mesh.elements->firstDomain; });

		clearIntervals();
		pushInterval(os.str().size());
		storeIntervals(name, os.str(), commitIntervals());
	}

	{ // CLUSTERS
		std::string filename = _name + ".CLUSTERS";
		std::string name = _path + filename;

		if (environment->MPIrank == 0) {
			(*_casefile) << "scalar per element:\tCLUSTERS\t\t" << filename << "\n";
			_casefile->flush();
		}

		std::stringstream os;

		if (environment->MPIrank == 0) {
			os << "CLUSTERS\n";
			os << "part\n";
			os << std::setw(10) << 1 << "\n";
		}

		eslocal cluster = _mesh.elements->gatherClustersDistribution()[environment->MPIrank];
		iterateElements(os, [&] (eslocal domain)->double { return _mesh.elements->clusters[domain] + cluster; });

		clearIntervals();
		pushInterval(os.str().size());
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
		os << "part\n";
		os << std::setw(10) << 1 << "\n";
		os << "coordinates\n";
	}

	os << std::showpos;
	for (size_t n = 0; n < _mesh.nodes->data.front()->gathredData->size(); ++n) {
		os << std::scientific << std::setprecision(5) << (*_mesh.nodes->data.front()->gathredData)[n] << "\n";
	}

	clearIntervals();
	pushInterval(os.str().size());
	storeIntervals(name, os.str(), commitIntervals());
}




;
