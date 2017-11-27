
#include "ensight.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/logging/logging.h"
#include "../../../basis/utilities/communication.h"
#include "../../../basis/utilities/utils.h"

#include "../../../config/ecf/environment.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"

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

		(*_casefile).flush();
	}
}

EnSight::~EnSight()
{
	if (_casefile != NULL) {
		delete _casefile;
	}
}

void EnSight::storeGeometry()
{
	std::string name = _path + _name + ".geo";
	std::vector<eslocal> ndistribution = _mesh.nodes->gatherUniqueNodeDistribution();
	std::vector<eslocal> edistribution = _mesh.elements->gatherElementsProcDistribution();
	eslocal loffset = 0, goffset = 0, lsize = 0, gsize = 0;

	std::vector<MPI_Aint> displacement;
	std::vector<int> lenghts;

	std::stringstream os;
	if (environment->MPIrank == 0) {
		os << "EnSight Gold geometry format\n";
		os << "----------------------------\n";

		os << "node id off\n";
		os << "element id off\n";

		os << "part\n";
		os << std::setw(10) << 1 << "\n";
		os << "MESH\n";

		os << "coordinates\n";
		os << std::setw(10) << ndistribution.back() << "\n";
	}

	auto pushInterval = [&] () {
		loffset = os.str().size() - lsize;
		lenghts.push_back(loffset);
		lsize = os.str().size();

		goffset = gsize;
		gsize += Communication::exscan(loffset);
		displacement.push_back(goffset + loffset);
	};

	auto storeNodes = [&] (std::function<double(const Point *p)> getCoordinate) {
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i) {
			if (_mesh.nodes->pintervals[i].sourceProcess == environment->MPIrank) {
				auto begin = _mesh.nodes->coordinates->datatarray().begin() + _mesh.nodes->pintervals[i].begin;
				auto end = _mesh.nodes->coordinates->datatarray().begin() + _mesh.nodes->pintervals[i].end;
				for (auto n = begin; n != end; ++n) {
					os << std::showpos << std::scientific << std::setprecision(5) << getCoordinate(n) << "\n";
				}
			}
		}
		pushInterval();
	};

	storeNodes([] (const Point *p) { return p->x; });
	storeNodes([] (const Point *p) { return p->y; });
	storeNodes([] (const Point *p) { return p->z; });

	if (environment->MPIrank == 0) {
		os << "hexa8\n";
		os << std::setw(10) << edistribution.back() << "\n";
	}
	auto node = _mesh.elements->nodes->datatarray().cbegin();
	for (size_t e = 0; e < _mesh.elements->size; e++) {
		for (size_t n = 0; n < 8; ++n, ++node) {
			auto it = std::lower_bound(_mesh.nodes->pintervals.begin(), _mesh.nodes->pintervals.end(), *node, [] (const ProcessInterval &interval, eslocal node) { return interval.end <= node; });
			os << std::setw(10) << it->globalOffset + *node - it->begin + 1;
		}
		os << "\n";
	}

	pushInterval();

	MPI_Datatype indexes;
	MPI_Type_create_hindexed(displacement.size(), lenghts.data(), displacement.data(), MPI_BYTE, &indexes);
	MPI_Type_commit(&indexes);

	MPI_File MPIfile;

	if (MPI_File_open(_storeCommunicator, name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MPIfile)) {
		ESINFO(ERROR) << "Unable to open file to store solution.";
	} else {
		MPI_File_set_view(MPIfile, 0, MPI_BYTE, indexes, "native", MPI_INFO_NULL);
		MPI_File_write_all(MPIfile, os.str().c_str(), os.str().size(), MPI_BYTE, MPI_STATUSES_IGNORE);
		MPI_File_close(&MPIfile);
	}
	MPI_Type_free(&indexes);
}

void EnSight::storeVariables()
{
	if (_mesh.nodes->data.size() != 1) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement store variables.";
	}

	std::string filename = _name + "." + _mesh.nodes->data.front()->names.front().substr(0, 4);
	std::string name = _path + filename;

	if (environment->MPIrank == 0) {
		(*_casefile) << "VARIABLE\n\n";
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

	for (size_t n = 0; n < _mesh.nodes->data.front()->gathredData->size(); ++n) {
		os << std::showpos << std::scientific << std::setprecision(5) << (*_mesh.nodes->data.front()->gathredData)[n] << "\n";
	}

	std::vector<MPI_Aint> displacement;
	std::vector<int> lenghts;

	eslocal offset = os.str().size();
	lenghts.push_back(offset);
	eslocal size = Communication::exscan(offset);
	displacement.push_back(offset);

	MPI_Datatype indexes;
	MPI_Type_create_hindexed(displacement.size(), lenghts.data(), displacement.data(), MPI_BYTE, &indexes);
	MPI_Type_commit(&indexes);

	MPI_File MPIfile;

	if (MPI_File_open(_storeCommunicator, name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MPIfile)) {
		ESINFO(ERROR) << "Unable to open file to store solution.";
	} else {
		MPI_File_set_view(MPIfile, 0, MPI_BYTE, indexes, "native", MPI_INFO_NULL);
		MPI_File_write_all(MPIfile, os.str().c_str(), os.str().size(), MPI_BYTE, MPI_STATUSES_IGNORE);
		MPI_File_close(&MPIfile);
	}
	MPI_Type_free(&indexes);

}




;
