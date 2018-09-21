
#ifndef SRC_INPUT_MPILOADER_MPILOADER_H_
#define SRC_INPUT_MPILOADER_MPILOADER_H_

#include <string>
#include <vector>

namespace espreso {

struct ParallelFile {
	const char *begin, *end;
	std::vector<char> data;
	std::vector<size_t> offsets;
};

struct MPILoader {

	static bool read(const std::string &file, ParallelFile &pfile, size_t align);
	static void align(ParallelFile &pfile, size_t lines);
};

}



#endif /* SRC_INPUT_MPILOADER_MPILOADER_H_ */
