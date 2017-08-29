
#include "mpi.h"

#include "../config/ecf/ecf.h"
#include "../config/reader/reader.h"
#include "../basis/logging/logging.h"

#include <iostream>
#include <fstream>

using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	ECFConfiguration ecf(&argc, &argv);
	ECFRedParameters redParameters = ECFReader::read(ecf, ECFReader::configurationFile, ecf.default_args, ecf.variables);
	std::ofstream os(ECFReader::configurationFile);
	ECFReader::store(ecf, os, true, false, redParameters);

	ESINFO(OVERVIEW) << "ECF is correct.";

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}



