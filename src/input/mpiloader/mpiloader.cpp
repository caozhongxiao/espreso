
#include "mpiloader.h"

#include "../../basis/containers/tarray.h"
#include "../../basis/logging/timeeval.h"
#include "../../basis/utilities/communication.h"

using namespace espreso;

bool MPILoader::read(const std::string &file, ParallelFile &pfile, size_t alignment)
{
	TimeEval timing("Read data from file");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("FILE OPEN");
	e1.start();

	MPI_File MPIfile;

	if (MPI_File_open(environment->MPICommunicator, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &MPIfile)) {
		return false;
	}

	e1.end();
	timing.addEvent(e1);

	TimeEvent e2("GET SIZE");
	e2.start();

	MPI_Offset size;
	MPI_File_get_size(MPIfile, &size);

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("STD VECTORS");
	e3.start();

	size_t block = 1;
	while (size / (1L << (block - 1)) > (1L << 31)) {
		++block;
	}
	block = 1L << block;

	std::vector<size_t> fdistribution = tarray<int>::distribute(environment->MPIsize, size / block + ((size % block) ? 1 : 0));
	std::vector<MPI_Aint> displacement = { (MPI_Aint)(block * fdistribution[environment->MPIrank]) };
	std::vector<int> length = { (int)(fdistribution[environment->MPIrank + 1] - fdistribution[environment->MPIrank]) };

	pfile.data.resize(block * length.front() + alignment);

	MPI_Datatype chunk;
	MPI_Datatype fDataDistribution;

	e3.end();
	timing.addEvent(e3);

	TimeEvent e4("COMMIT DATATYPES");
	e4.start();

	MPI_Type_contiguous(block, MPI_BYTE, &chunk);
	MPI_Type_commit(&chunk);
	MPI_Type_create_hindexed(1, length.data(), displacement.data(), chunk, &fDataDistribution);
	MPI_Type_commit(&fDataDistribution);

	e4.end();
	timing.addEvent(e4);

	TimeEvent e5("SET VIEW");
	e5.start();

	MPI_File_set_view(MPIfile, 0, chunk, fDataDistribution, "native", MPI_INFO_NULL);

	e5.end();
	timing.addEvent(e5);

	TimeEvent e6("READ ALL");
	e6.start();

	MPI_File_read_all(MPIfile, pfile.data.data(), length.front(), chunk, MPI_STATUS_IGNORE);

	e6.end();
	timing.addEvent(e6);

	MPI_File_close(&MPIfile);
	MPI_Type_free(&chunk);

	pfile.begin = pfile.data.data();
	pfile.end = pfile.data.data() + block * length.front();
	if (environment->MPIrank + 1 == environment->MPIsize) {
		pfile.end -= block * fdistribution.back() - size;
	}
	pfile.offsets = Communication::getDistribution<size_t>(pfile.end - pfile.begin);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();

	return true;
}

void MPILoader::align(ParallelFile &pfile, size_t lines)
{
	const char *_current = pfile.data.data();
	if (environment->MPIsize > 1) {
		const char *exchangeStart = _current;
		if (environment->MPIrank) {
			while (*_current++ != '\n');
			exchangeStart = _current;
			for (int i = 1; i < lines; i++) {
				while (*exchangeStart++ != '\n');
			}
		}

		std::vector<std::vector<char> > sBuffer(1, std::vector<char>(const_cast<const char*>(pfile.data.data()), exchangeStart)), rBuffer;

		std::vector<int> neighs;
		if (environment->MPIrank) {
			neighs.push_back(environment->MPIrank - 1);
		}
		if (environment->MPIrank + 1 < environment->MPIsize) {
			neighs.push_back(environment->MPIrank + 1);
		}
		rBuffer.resize(neighs.size());
		Communication::receiveUpperUnknownSize(sBuffer, rBuffer, neighs);

		if (rBuffer.back().size() > pfile.data.size() - (pfile.end - pfile.begin)) {
			pfile.data.resize((pfile.end - pfile.begin) + rBuffer.back().size());
		}
		memcpy(pfile.data.data() + (pfile.end - pfile.begin), rBuffer.back().data(), rBuffer.back().size());
		char* firstLineEnd = rBuffer.back().data();
		if (rBuffer.back().size()) {
			while (*firstLineEnd++ != '\n');
		}
		pfile.end += firstLineEnd - rBuffer.back().data();
	}
	pfile.begin = _current;
	pfile.offsets = Communication::getDistribution<size_t>(pfile.end - pfile.begin);
}

