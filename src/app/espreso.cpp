
#include "esinfo/eslog.hpp"
#include "esinfo/timeinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/utilities/communication.h"

#include "config/ecf/root.h"
#include "config/reader/reader.h"
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
	info::env::set();
	info::mpi::set();

	ECFReader ecf(argc, argv);

	bool divided = info::mpi::divide(ecf.meshDuplication); // we cannot print error before initialization of the logger
	eslog::init(ecf.ecf.c_str(), ecf.outputPath.c_str(), ecf.meshDuplication);
	if (!divided) {
		eslog::globalerror("Cannot set MESH DUPLICATION: the number of MPI processes is not divisible by %d\n", ecf.meshDuplication);
	}

	eslog::start("ESPRESO: STARTED", "ESPRESO");
	eslog::param("MPI", info::mpi::size);
	eslog::param("OMP/MPI", info::env::OMP_NUM_THREADS);
	if (ecf.meshDuplication > 1) {
		eslog::param("MESH DUPLICATION", ecf.meshDuplication);
	}
	eslog::ln();

	Mesh::init();

	info::ecf = new ECFRoot(&argc, &argv);
	info::mesh = new Mesh();

	MPITools::init();

	info::mesh->store = ResultStore::createAsynchronizedStore(*info::mesh);

	eslog::checkpoint("ESPRESO: CONFIGURED");
	eslog::param("ecf", info::ecf->ecffile.c_str());
	eslog::ln();

	if (ResultStore::isComputeNode()) {
		Input::load(*info::ecf, *info::mesh);
		eslog::checkpoint("ESPRESO: MESH LOADED");
		eslog::param("database", Input::inputFile(*info::ecf));
		eslog::ln();
		if (info::mpi::irank == 0) {
			info::mesh->preprocess();
			eslog::checkpointln("ESPRESO: MESH PREPROCESSED");
			info::mesh->printMeshStatistics();
		}
		if (info::mpi::isize > 1) {
			info::mesh->duply();
			eslog::checkpointln("ESPRESO: MESH DUPLICATED");
		}

		if (Input::convertDatabase(*info::ecf)) {
			info::mesh->storeMesh();
			eslog::endln("ESPRESO: MESH STORED");
		} else {
			LoadStepIterator steps;
			while (steps.next()) {
				eslog::checkpoint("ESPRESO: SOLVED");
				eslog::param("LOADSTEP", time::step);
				eslog::ln();
			}
			eslog::end("ESPRESO: SOLVED");
			eslog::param("LOADSTEP", time::step);
			eslog::ln();
		}
	}

	eslog::finish();
	ResultStore::destroyAsynchronizedStore();
	Mesh::destroy();
	MPITools::destroy();
	info::mpi::finish();

	delete info::ecf;
	delete info::mesh;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
