
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/utilities/communication.h"
#include "basis/logging/logging.h"

#include "config/ecf/root.h"
#include "physics/loadstepiterator.h"

#include "mesh/mesh.h"
#include "input/input.h"
#include "output/result/resultstore.h"

using namespace espreso;

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	info::system::setSignals();
	Logging::init();
	Mesh::init();
	MPITools::init();
	info::env::set();
	info::mpi::set();

	info::ecf = new ECFRoot(&argc, &argv);
	info::mesh = new Mesh();
	info::mesh->store = ResultStore::createAsynchronizedStore(*info::mesh);

	ESINFO(OVERVIEW) <<
			"Starting ESPRESO, " <<
			"MPI: " << info::mpi::size << ", "
			"OMP/MPI: " << info::env::OMP_NUM_THREADS;

	if (ResultStore::isComputeNode()) {
		if (Input::load(*info::ecf, *info::mesh)) {
			LoadStepIterator steps;
			while (steps.next());
		} else {
			info::mesh->storeMesh();
		}
	}

	ResultStore::destroyAsynchronizedStore();
	Mesh::destroy();
	MPITools::destroy();
	delete info::ecf;
	delete info::mesh;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
