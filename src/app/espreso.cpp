
#include "espreso.h"

#include "../basis/logging/logging.h"

#include "../config/ecf/root.h"
#include "../globals/env.h"
#include "../globals/system.h"

//#include "../physics/loadstepiterator.h"

#include "../mesh/mesh.h"
#include "../input/input.h"
#include "../output/result/resultstore.h"

int espreso::run(int *argc, char ***argv)
{
	int provided;
	MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);

	system::setSignals();
	env::setMPI();

	ECFRoot ecf(argc, argv);
	Mesh mesh(ecf, NULL);
	mesh.store = ResultStore::createAsynchronizedStore(mesh, ecf.output);


	std::string processes, threads;
	if (mesh.store->storeProcesses) {
		processes = std::to_string(mesh.store->computeProcesses) + " + " + std::to_string(mesh.store->storeProcesses);
	} else {
		processes = std::to_string(mesh.store->computeProcesses);
	}
	if (mesh.store->storeThreads) {
		threads = std::to_string(environment->OMP_NUM_THREADS) + " + " + std::to_string(mesh.store->storeThreads);
	} else {
		threads = std::to_string(environment->OMP_NUM_THREADS);
	}

	ESINFO(OVERVIEW) << "Run ESPRESO SOLVER using " << processes << " MPI and " << threads << " threads.";

	auto computeSolution = [&] () {
		switch (ecf.input) {
		case INPUT_FORMAT::WORKBENCH:
			return !ecf.workbench.convert_database;
		case INPUT_FORMAT::OPENFOAM:
			return !ecf.openfoam.convert_database;
		default:
			return true;
		}
	};

	if (ResultStore::isComputeNode()) {

		Input::load(ecf, mesh);

//		LoadStepIterator steps(mesh, ecf.heat_transfer_2d);

//		for (size_t step = 0; step < _loader->loadSteps(); step++) {
//			_loadSteps.push_back(_loader->getLoadStepSolver(step, _mesh, _store));
//		}
//
//		_loader->preprocessMesh();
//		mesh.initNodeData();
//
//		mesh.preprocessing->finishPreprocessing();
//
//		_store->updateMesh();
//
//
//		Factory factory(ecf, mesh, *solutionStore);
//		if (computeSolution()) {
//			factory.solve();
//		} else {
//			// ESPRESOBinaryFormat::store(mesh, configuration);
//		}
	}

	ResultStore::destroyAsynchronizedStore();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}




