
#ifndef SRC_INPUT_MPILOADER_MPILOADER_H_
#define SRC_INPUT_MPILOADER_MPILOADER_H_

#include <string>
#include <vector>

#include "../../basis/utilities/communication.h"

namespace espreso {

struct LoaderConfiguration;

struct ParallelFile {
	const char *begin, *end;
	std::vector<char> data;
	std::vector<size_t> offsets;
};

struct MPILoader {

	static bool open(MPITools::MPIGroup &group, MPI_File &MPIfile, const std::string &file);
	static void read(MPITools::MPIGroup &group, MPI_File &MPIfile, ParallelFile &pfile, size_t alignment);
	static void scatter(MPITools::MPIGroup &group, ParallelFile &pfile, size_t alignment);
	static void align(MPITools::MPIGroup &group, ParallelFile &pfile, size_t lines);
};

}



#endif /* SRC_INPUT_MPILOADER_MPILOADER_H_ */
