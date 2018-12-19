
#include "../basis/logging/logging.h"

#include "../config/ecf/root.h"
#include "../globals/run.h"
#include "../globals/env.h"
#include "../globals/system.h"

#include "../physics/loadstepiterator.h"

#include "../mesh/mesh.h"
#include "../input/input.h"
#include "../output/result/resultstore.h"

using namespace espreso;

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	system::setSignals();
	env::setMPI();

	run::ecf = new ECFRoot(&argc, &argv);
	run::mesh = new Mesh();
	run::mesh->store = ResultStore::createAsynchronizedStore(*run::mesh, run::ecf->output);

	ESINFO(OVERVIEW) <<
			"Starting ESPRESO, " <<
			"MPI: " << environment->MPIsize << ", "
			"OMP/MPI: " << environment->OMP_NUM_THREADS;

	if (ResultStore::isComputeNode()) {

		Input::load(*run::ecf, *run::mesh);

		LoadStepIterator steps;
		while (steps.next());
	}

	ResultStore::destroyAsynchronizedStore();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
