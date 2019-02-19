

#include "ensight.h"

#include "esinfo/timeinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"

#include "mesh/elements/element.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <algorithm>
#include <functional>
#include <fstream>


using namespace espreso;

EnSight::EnSight(const std::string &name, const Mesh &mesh)
: CollectedVisualization(mesh), _path(std::string(eslog::path()) + "/"), _name(name), _variableCounter(0)
{
	_caseheader << "#\n";
	_caseheader << "# ESPRESO solution\n";
	_caseheader << "#\n";

	_caseheader << "\nFORMAT\n";
	_caseheader << "type: \tensight gold\n\n";

	_casegeometry << "GEOMETRY\n\n";

	_casevariables << "VARIABLE\n\n";
}

EnSight::~EnSight()
{

}

void EnSight::storecasefile()
{
	if (info::mpi::rank == 0) {
		std::ofstream os(_path + _name + ".case");
		os << _caseheader.str() << "\n";
		os << _casegeometry.str() << "\n";
		os << _casevariables.str() << "\n";

		if (_variableCounter) {
			os << "TIME\n\n";
			os << "time set:              1\n";
			os << "number of steps:       " << _variableCounter << "\n";
			os << "filename start number: 1\n";
			os << "filename increment:    1\n";
			os << "time values:           " << _casetime.str() << "\n";
		}
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
		eslog::error("ESPRESO internal error: unknown element code.\n");
		return "";
	}
}

void EnSight::updateMesh()
{
	_casegeometry << "model:\t" << _directory << _name << ".geo\n\n";

	std::string name = _path + _directory + _name + ".geo";
	int part = 1;

	std::stringstream os;
	os << std::showpos << std::scientific << std::setprecision(5);

	if (info::mpi::rank == 0) {
		_writer.storeFormat(os);
		_writer.storeDescriptionLine(os, "EnSight Gold geometry format");
		_writer.storeDescriptionLine(os, "----------------------------");

		_writer.storeDescriptionLine(os, "node id off");
		_writer.storeDescriptionLine(os, "element id off");
	}

	auto storePartHeader = [&] (const std::string &name, esint nodes) {
		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, "part");
			_writer.storeInt(os, part++);
			_writer.storeDescriptionLine(os, name);

			_writer.storeDescriptionLine(os, "coordinates");
			_writer.storeInt(os, nodes);
		}
	};

	auto storeRegionNodes = [&] (const std::vector<ProcessInterval> &intervals, const serializededata<esint, esint> *nodes, std::function<double(const Point *p)> getCoordinate) {
		for (size_t i = 0; i < intervals.size(); i++) {
			if (intervals[i].sourceProcess == info::mpi::rank) {
				for (auto n = nodes->datatarray().cbegin() + intervals[i].begin; n != nodes->datatarray().cbegin() + intervals[i].end; ++n) {
					_writer.storeFloat(os, getCoordinate(_mesh.nodes->coordinates->datatarray().cbegin() + *n));
				}
			}
		}
		pushInterval(os.str().size());
	};

	auto storeNodeRegionIndex = [&] (esint n, const std::vector<ProcessInterval> &nintervals, const tarray<esint> &nodes) {
		auto iit = std::lower_bound(_mesh.nodes->pintervals.begin(), _mesh.nodes->pintervals.end(), n, [] (const ProcessInterval &interval, esint node) { return interval.end <= node; });
		size_t iindex = iit - _mesh.nodes->pintervals.begin();
		esint offset = std::lower_bound(nodes.begin() + nintervals[iindex].begin, nodes.begin() + nintervals[iindex].end, n) - nodes.begin();
		return nintervals[iindex].globalOffset + offset - nintervals[iindex].begin;
	};

	auto storeElements = [&] (
			const std::vector<esint> &ecounters,
			const tarray<esint> &elements, const std::vector<ElementsInterval> &eintervals,
			const tarray<esint> &nodes, const std::vector<ProcessInterval> &nintervals) {

		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (ecounters[etype]) {
				if (info::mpi::rank == 0) {
					_writer.storeDescriptionLine(os, codetotype(etype));
					_writer.storeInt(os, ecounters[etype]);
				}

				for (size_t i = 0; i < eintervals.size(); i++) {
					if (eintervals[i].code == etype) {
						auto enodes = _mesh.elements->procNodes->cbegin() + eintervals[i].begin;
						esint prev = eintervals[i].begin;
						for (esint e = eintervals[i].begin; e < eintervals[i].end; prev = elements[e++]) {
							enodes += elements[e] - prev;
							for (auto n = enodes->begin(); n != enodes->end(); ++n) {
								_writer.eIndex(os, storeNodeRegionIndex(*n, nintervals, nodes) + 1);
							}
							_writer.eEnd(os);
						}
					}
				}
				pushInterval(os.str().size());
			}
		}
	};

	auto storeBoundaryElements = [&] (const BoundaryRegionStore *region) {
		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (region->ecounters[etype]) {
				if (info::mpi::rank == 0) {
					_writer.storeDescriptionLine(os, codetotype(etype));
					_writer.storeInt(os, region->ecounters[etype]);
				}

				for (size_t i = 0; i < region->eintervals.size(); i++) {
					if (region->eintervals[i].code == etype) {
						auto enodes = region->procNodes->cbegin() + region->eintervals[i].begin;
						for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e, ++enodes) {
							for (auto n = enodes->begin(); n != enodes->end(); ++n) {
								_writer.eIndex(os, storeNodeRegionIndex(*n, region->nintervals, region->nodes->datatarray()) + 1);
							}
							_writer.eEnd(os);
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

		storeElements(region->ecounters, region->elements->datatarray(), region->eintervals, region->nodes->datatarray(), region->nintervals);
	};

	auto storeBRegion = [&] (const BoundaryRegionStore *region) {
		storeRegionNodes(region->nintervals, region->nodes, [] (const Point *p) { return p->x; });
		storeRegionNodes(region->nintervals, region->nodes, [] (const Point *p) { return p->y; });
		storeRegionNodes(region->nintervals, region->nodes, [] (const Point *p) { return p->z; });

		if (region->dimension) {
			storeBoundaryElements(region);
		} else {
			if (info::mpi::rank == 0) {
				_writer.storeDescriptionLine(os, codetotype(static_cast<int>(Element::CODE::POINT1)));
				_writer.storeInt(os, region->uniqueTotalSize);
			}
			for (esint i = 0; i < region->uniqueSize; ++i) {
				_writer.storeInt(os, region->uniqueOffset + i + 1);
			}
			pushInterval(os.str().size());
		}
	};

	for (size_t r = 1; r < _mesh.elementsRegions.size(); r++) {
		storePartHeader(_mesh.elementsRegions[r]->name, _mesh.elementsRegions[r]->uniqueTotalSize);
		storeERegion(_mesh.elementsRegions[r]);
	}

	for (size_t r = 1; r < _mesh.boundaryRegions.size(); r++) {
		storePartHeader(_mesh.boundaryRegions[r]->name, _mesh.boundaryRegions[r]->uniqueTotalSize);
		storeBRegion(_mesh.boundaryRegions[r]);
	}

	storeIntervals(name, os.str(), commitIntervals());

	storecasefile();
}

void EnSight::setvariables()
{
	auto tabs = [] (size_t size) {
		return size < 8 ? 2 : 1;
	};

	auto format = [] (size_t size) {
		if (size == 1) return "scalar";
		if (size == 3 || size == 2) return "vector";
		return "";
	};

	for (size_t i = 0; i < _mesh.nodes->data.size(); i++) {
		if (_mesh.nodes->data[i]->names.size()) {
			std::string filename = _directory + _mesh.nodes->data[i]->names.front();
			_casevariables << format(_mesh.nodes->data[i]->dimension);
			_casevariables << " per node:\t1 ";
			_casevariables << _mesh.nodes->data[i]->names.front();
			_casevariables << std::string(tabs(_mesh.nodes->data[i]->names.front().size()), '\t');
			_casevariables << filename << ".****\n";
		}
	}
	for (size_t i = 0; i < _mesh.elements->data.size(); i++) {
		if (_mesh.elements->data[i]->names.size()) {
			std::string filename = _directory + _mesh.elements->data[i]->names.front();
			_casevariables << format(_mesh.elements->data[i]->dimension);
			_casevariables << " per element:\t1 ";
			_casevariables << _mesh.elements->data[i]->names.front();
			_casevariables << std::string(tabs(_mesh.elements->data[i]->names.front().size()), '\t');
			_casevariables << filename << ".****\n";
		}
	}
}

void EnSight::storeDecomposition()
{
	_casevariables << "scalar per element:\tDOMAINS\t\t" << _directory << "DOMAINS" << "\n";
	_casevariables << "scalar per element:\tCLUSTERS\t" << _directory << "CLUSTERS" << "\n";
	_casevariables << "scalar per element:\tMPI\t" << _directory << "MPI" << "\n";

	auto iterateElements = [&] (std::stringstream &os, const std::vector<ElementsInterval> &intervals, const std::vector<esint> &ecounters, std::function<double(esint domain)> fnc) {
		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (ecounters[etype]) {
				if (info::mpi::rank == 0) {
					_writer.storeDescriptionLine(os, codetotype(etype));
				}

				for (size_t i = 0; i < intervals.size(); i++) {
					if (intervals[i].code == etype) {
						for (esint e = intervals[i].begin; e < intervals[i].end; ++e) {
							_writer.storeFloat(os, fnc(intervals[i].domain));
						}
					}
				}

				pushInterval(os.str().size());
			}
		}
	};

	esint part = 1;

	auto storePartHeader = [&] (std::stringstream &os) {
		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, "part");
			_writer.storeInt(os, part++);
		}
	};


	{ // DOMAINS
		std::string filename = _directory + "DOMAINS";
		std::string name = _path + filename;

		std::stringstream os;
		os << std::showpos << std::scientific << std::setprecision(5);

		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, "DOMAINS");
		}

		part = 1;
		clearIntervals();
		for (size_t r = 1; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader(os);
			iterateElements(os, _mesh.elementsRegions[r]->eintervals, _mesh.elementsRegions[r]->ecounters, [&] (esint domain)->double { return domain; });
		}
		storeIntervals(name, os.str(), commitIntervals());
	}

	{ // CLUSTERS
		std::string filename = _directory + "CLUSTERS";
		std::string name = _path + filename;

		std::stringstream os;
		os << std::showpos << std::scientific << std::setprecision(5);

		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, "CLUSTERS");
		}

		part = 1;
		clearIntervals();
		esint cluster = _mesh.elements->gatherClustersDistribution()[info::mpi::rank];
		for (size_t r = 1; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader(os);
			iterateElements(os, _mesh.elementsRegions[r]->eintervals, _mesh.elementsRegions[r]->ecounters, [&] (esint domain)->double { return _mesh.elements->clusters[domain - _mesh.elements->firstDomain] + cluster; });
		}

		storeIntervals(name, os.str(), commitIntervals());
	}

	{ // MPI
		std::string filename = _directory + "MPI";
		std::string name = _path + filename;

		std::stringstream os;
		os << std::showpos << std::scientific << std::setprecision(5);

		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, "MPI");
		}

		part = 1;
		clearIntervals();
		for (size_t r = 1; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader(os);
			iterateElements(os, _mesh.elementsRegions[r]->eintervals, _mesh.elementsRegions[r]->ecounters, [&] (esint domain)->double { return info::mpi::rank; });
		}

		storeIntervals(name, os.str(), commitIntervals());
	}

	storecasefile();
}

void EnSight::updateSolution()
{
	if (!Visualization::storeStep()) {
		return;
	}

	if (_variableCounter == 0) {
		setvariables();
	}

	_casetime << time::current;
	if ((++_variableCounter) % 10 == 0) {
		_casetime << "\n                       ";
	} else {
		_casetime << " ";
	}

	for (size_t di = 0; di < _mesh.nodes->data.size(); di++) {
		if (_mesh.nodes->data[di]->names.size() == 0) {
			continue;
		}
//		Communication::serialize([&] () {
//			std::cout << _mesh.nodes->data[di]->data;
//		});
		esint size = _mesh.nodes->data[di]->dimension;

		std::string filename = _directory + _mesh.nodes->data[di]->names.front();
		std::stringstream name;
		name << _path + filename + "." << std::setw(4) << std::setfill('0') << _variableCounter;

		std::stringstream os;
		os << std::showpos << std::scientific << std::setprecision(5);

		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, _mesh.nodes->data[di]->names.front());
		}

		esint part = 1;

		auto storePartHeader = [&] () {
			if (info::mpi::rank == 0) {
				_writer.storeDescriptionLine(os, "part");
				_writer.storeInt(os, part++);
				_writer.storeDescriptionLine(os, "coordinates");
			}
		};

		auto iterateNodes = [&] (const std::vector<ProcessInterval> &intervals, const tarray<esint> &nodes) {
			for (esint s = 0; s < size; s++) {
				for (size_t i = 0; i < intervals.size(); i++) {
					if (intervals[i].sourceProcess == info::mpi::rank) {
						esint offset = _mesh.nodes->pintervals[i].globalOffset - _mesh.nodes->uniqueOffset + (_mesh.nodes->size - _mesh.nodes->uniqueSize);
						for (esint n = intervals[i].begin; n < intervals[i].end; ++n) {
							esint index = offset + nodes[n] - _mesh.nodes->pintervals[i].begin;
							_writer.storeFloat(os, _mesh.nodes->data[di]->data[size * index + s]);
						}
					}
				}
				pushInterval(os.str().size());
			}
			if (size == 2) {
				for (size_t i = 0; i < intervals.size(); i++) {
					if (intervals[i].sourceProcess == info::mpi::rank) {
						for (esint n = intervals[i].begin; n < intervals[i].end; ++n) {
							_writer.storeFloat(os, .0);
						}
					}
				}
				pushInterval(os.str().size());
			}
		};

		clearIntervals();
		for (size_t r = 1; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader();
			iterateNodes(_mesh.elementsRegions[r]->nintervals, _mesh.elementsRegions[r]->nodes->datatarray());
		}

		for (size_t r = 1; r < _mesh.boundaryRegions.size(); r++) {
			storePartHeader();
			iterateNodes(_mesh.boundaryRegions[r]->nintervals, _mesh.boundaryRegions[r]->nodes->datatarray());
		}

		storeIntervals(name.str(), os.str(), commitIntervals());
	}

	for (size_t di = 0; di < _mesh.elements->data.size(); di++) {
		if (_mesh.elements->data[di]->names.size() == 0) {
			continue;
		}
		esint size = _mesh.elements->data[di]->dimension;

		std::string filename = _directory + _mesh.elements->data[di]->names.front();
		std::stringstream name;
		name << _path + filename + "." << std::setw(4) << std::setfill('0') << _variableCounter;

		std::stringstream os;
		os << std::showpos << std::scientific << std::setprecision(5);

		if (info::mpi::rank == 0) {
			_writer.storeDescriptionLine(os, _mesh.elements->data[di]->names.front());
		}

		esint part = 1;

		auto storePartHeader = [&] () {
			if (info::mpi::rank == 0) {
				_writer.storeDescriptionLine(os, "part");
				_writer.storeInt(os, part++);
			}
		};

		auto iterateElements = [&] (std::stringstream &os, int etype, const tarray<esint> &elements, const std::vector<ElementsInterval> &eintervals) {
			for (esint s = 0; s < size; s++) {
				for (size_t i = 0; i < eintervals.size(); i++) {
					if (eintervals[i].code == etype) {
						for (esint e = eintervals[i].begin; e < eintervals[i].end; ++e) {
							_writer.storeFloat(os, _mesh.elements->data[di]->data[size * elements[e] + s]);
						}
					}
				}
				pushInterval(os.str().size());
			}
			if (size == 2) {
				for (size_t i = 0; i < eintervals.size(); i++) {
					if (eintervals[i].code == etype) {
						for (esint e = eintervals[i].begin; e < eintervals[i].end; ++e) {
							_writer.storeFloat(os, .0);
						}
					}
				}
				pushInterval(os.str().size());
			}
		};

		clearIntervals();
		for (size_t r = 1; r < _mesh.elementsRegions.size(); r++) {
			storePartHeader();
			for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
				if (_mesh.elementsRegions[r]->ecounters[etype]) {
					if (info::mpi::rank == 0) {
						_writer.storeDescriptionLine(os, codetotype(etype));
					}
					iterateElements(os, etype, _mesh.elementsRegions[r]->elements->datatarray(), _mesh.elementsRegions[r]->eintervals);
				}
			}
		}

		storeIntervals(name.str(), os.str(), commitIntervals());
	}

	storecasefile();
}



