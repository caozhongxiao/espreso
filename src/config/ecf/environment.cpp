
#include "environment.h"
#include "../configuration.hpp"
#include "../../basis/utilities/utils.h"

namespace espreso {

Environment *environment = NULL;

Environment::Environment(): executable("espreso")
{
	int initialized;
	MPI_Initialized(&initialized);

	MPICommunicator = MPI_COMM_WORLD;
	if (initialized) {
		MPI_Comm_rank(MPICommunicator, &MPIrank);
		MPI_Comm_size(MPICommunicator, &MPIsize);
	}

	char *var = getenv("_");
	if (var != NULL) {
		executable = Esutils::getEnv<std::string>("_");
	} else {
		executable = "espreso";
	}

	MKL_NUM_THREADS    = Esutils::getEnv<size_t>("MKL_NUM_THREADS");
	OMP_NUM_THREADS    = Esutils::getEnv<size_t>("OMP_NUM_THREADS");
	SOLVER_NUM_THREADS = Esutils::getEnv<size_t>("SOLVER_NUM_THREADS");
	PAR_NUM_THREADS    = Esutils::getEnv<size_t>("PAR_NUM_THREADS");
	CILK_NWORKERS      = Esutils::getEnv<size_t>("CILK_NWORKERS");

	MPI_Bcast(&Logging::time, sizeof(time_t), MPI_BYTE, 0, MPICommunicator);

	log_dir = "debug";
	REGISTER(log_dir, ECFMetaData()
			.setdescription({ "A name of logging directory" })
			.setdatatype({ ECFDataType::STRING }));

	verbose_level = 1;
	measure_level = testing_level = 0;
	remove_old_results = print_matrices = false;

	REGISTER(verbose_level, ECFMetaData()
			.setdescription({ "Verbose level [0-3]." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

//	REGISTER(testing_level, ECFMetaData()
//			.setdescription({ "Testing level [0-3]." })
//			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	REGISTER(measure_level, ECFMetaData()
			.setdescription({ "Measure level [0-3]." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	REGISTER(print_matrices, ECFMetaData()
			.setdescription({ "Print assembler matrices for debugging." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(remove_old_results, ECFMetaData()
			.setdescription({ "Keep results of only the last ESPRESO run." })
			.setdatatype({ ECFDataType::BOOL }));

	if (environment == NULL) {
		environment = this;
	}
}

}



