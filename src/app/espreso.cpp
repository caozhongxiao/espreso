
#include <signal.h>
#include <csignal>

#include "mpi.h"

#include "factory/factory.h"
#include "../config/ecf/ecf.h"
#include "../config/reader/reader.h"
#include "../basis/logging/logging.h"

using namespace espreso;

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGSEGV:
		ESINFO(ERROR) << "Invalid memory reference";
		break;
	case SIGFPE:
		ESINFO(ERROR) << "Erroneous arithmetic operation";
		break;
	}
}


int main(int argc, char **argv)
{
	std::signal(SIGFPE, signalHandler);
	std::signal(SIGSEGV, signalHandler);

	MPI_Init(&argc, &argv);

	ECFConfiguration ecf(&argc, &argv);
	exit(0);

	ESINFO(OVERVIEW) << "Run ESPRESO on " << environment->MPIsize << " process(es).";

	Factory factory(ecf);

	factory.solve();
	factory.finalize();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}


