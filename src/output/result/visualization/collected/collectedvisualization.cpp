
#include "collectedvisualization.h"

#include "basis/utilities/communication.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

CollectedVisualization::CollectedVisualization(const Mesh &mesh)
: Visualization(mesh)
{
	MPI_Comm_split(info::mpi::comm, 0, info::mpi::rank, &_storeCommunicator);
	clearIntervals();
}

CollectedVisualization::~CollectedVisualization()
{
	MPI_Comm_free(&_storeCommunicator);
	for (size_t i = 0; i < _datatypes.size(); i++) {
		MPI_Type_free(_datatypes[i]);
		delete _datatypes[i];
	}
}

void CollectedVisualization::clearIntervals()
{
	_loffset = _goffset = _lsize = _gsize = 0;
	_displacement.clear();
	_lenghts.clear();
}

void CollectedVisualization::pushInterval(esint size)
{
	_loffset = size - _lsize;
	_lenghts.push_back(_loffset);
	_lsize = size;

	_goffset = _gsize;
	_gsize += Communication::exscan(_loffset);
	_displacement.push_back(_goffset + _loffset);
}

MPI_Datatype* CollectedVisualization::commitIntervals()
{
	_datatypes.push_back(new MPI_Datatype());
	MPI_Type_create_hindexed(_displacement.size(), _lenghts.data(), _displacement.data(), MPI_BYTE, _datatypes.back());
	MPI_Type_commit(_datatypes.back());
	return _datatypes.back();
}

void CollectedVisualization::storeIntervals(const std::string &name, const std::string &data, MPI_Datatype* datatype)
{
	MPI_File MPIfile;

	if (MPI_File_open(_storeCommunicator, name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MPIfile)) {
		eslog::error("MPI cannot create file '%s'\n", name.c_str());
	} else {
		MPI_File_set_view(MPIfile, 0, MPI_BYTE, *datatype, "native", MPI_INFO_NULL);
		MPI_File_write_all(MPIfile, data.c_str(), data.size(), MPI_BYTE, MPI_STATUS_IGNORE);
		MPI_File_close(&MPIfile);
	}
}
