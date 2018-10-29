
#include "mpiloader.h"

#include "../../basis/containers/tarray.h"

using namespace espreso;

bool MPILoader::open(MPIGroup &group, MPI_File &MPIfile, const std::string &file)
{
//	MPI_Info info;
//	MPI_Info_create(&info);
//	MPI_Info_set(info, "access_style", "read_once");
//	MPI_Info_set(info, "collective_buffering", "true");
//	MPI_Info_set(info, "romio_cb_read", "enable");

	if (MPI_File_open(group.communicator, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &MPIfile)) {
		return false;
	}

//	MPI_Info_free(&info);

	return true;
}

void MPILoader::read(MPIGroup &group, MPI_File &MPIfile, ParallelFile &pfile)
{
	MPI_Offset size;
	MPI_File_get_size(MPIfile, &size);

	size_t block = 1;
	while (size / (1L << (block - 1)) > (1L << 31)) {
		++block;
	}
	block = 1L << block;

	std::vector<int> fdistribution = tarray<int>::distribute(group.size, size / block + ((size % block) ? 1 : 0));
	std::vector<MPI_Aint> displacement = { (MPI_Aint)(block * fdistribution[group.rank]) };
	std::vector<int> length = { (int)(fdistribution[group.rank + 1] - fdistribution[group.rank]) };

	pfile.data.resize(block * length.front() + pfile.align);

	MPI_Datatype chunk, fDataDistribution;

	MPI_Type_contiguous(block, MPI_BYTE, &chunk);
	MPI_Type_commit(&chunk);
	MPI_Type_create_hindexed(1, length.data(), displacement.data(), chunk, &fDataDistribution);
	MPI_Type_commit(&fDataDistribution);

//	MPI_Info info;
//	MPI_Info_create(&info);
//	MPI_Info_set(info, "access_style", "read_once");
//	MPI_Info_set(info, "collective_buffering", "true");
//	MPI_Info_set(info, "romio_cb_read", "enable");

	MPI_File_set_view(MPIfile, 0, chunk, fDataDistribution, "native", MPI_INFO_NULL);
	MPI_File_read_all(MPIfile, pfile.data.data(), length.front(), chunk, MPI_STATUS_IGNORE);

//	MPI_Info_free(&info);
	MPI_File_close(&MPIfile);
	MPI_Type_free(&chunk);

	pfile.begin = pfile.data.data();
	pfile.end = pfile.data.data() + block * length.front();
	if (group.rank + 1 == group.size) {
		pfile.end -= block * fdistribution.back() - size;
	}
	pfile.offsets = Communication::getDistribution<size_t>(pfile.end - pfile.begin, group);
}

void MPILoader::scatter(MPIGroup &group, ParallelFile &pfile)
{
	if (group.size == 1) {
		return;
	}

	size_t datasize = pfile.data.size() - pfile.align;
	MPI_Bcast(&datasize, sizeof(size_t), MPI_BYTE, 0, group.communicator);
	size_t chunk = std::ceil((double)datasize / group.size);

	if (group.rank) {
		pfile.data.resize(chunk + pfile.align);
	} else if (pfile.data.size() < group.size * chunk) {
		pfile.data.resize(group.size * chunk);
	}

	std::vector<char> data(chunk + pfile.align);
	MPI_Scatter(pfile.data.data(), chunk, MPI_BYTE, data.data(), chunk, MPI_BYTE, 0, group.communicator);
	pfile.data.swap(data);

	pfile.begin = pfile.data.data();
	pfile.end = pfile.data.data() + chunk;
	if (group.rank + 1 == group.size) {
		pfile.end -= chunk * group.size - datasize;
	}

	pfile.offsets = Communication::getDistribution<size_t>(pfile.end - pfile.begin);
}

void MPILoader::bcast(MPIGroup &group, ParallelFile &pfile)
{
	size_t size = pfile.data.size();
	MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, 0, group.communicator);
	pfile.data.resize(size);
	MPI_Bcast(pfile.data.data(), pfile.data.size(), MPI_BYTE, 0, group.communicator);
	pfile.begin = pfile.data.data();
	pfile.end = pfile.data.data() + pfile.data.size();
}

void MPILoader::align(MPIGroup &group, ParallelFile &pfile, size_t lines)
{
	if (pfile.offsets.back() == 0) {
		return;
	}
	const char *_current = pfile.data.data();
	if (group.size > 1) {
		const char *exchangeStart = _current;
		if (group.rank) {
			while (*_current++ != '\n');
			exchangeStart = _current;
			for (size_t i = 1; i < lines; i++) {
				while (*exchangeStart++ != '\n');
			}
		}

		std::vector<std::vector<char> > sBuffer(1, std::vector<char>(const_cast<const char*>(pfile.data.data()), exchangeStart)), rBuffer;

		std::vector<int> neighs;
		if (group.rank) {
			neighs.push_back(group.rank - 1);
		}
		if (group.rank + 1 < group.size) {
			neighs.push_back(group.rank + 1);
		}
		rBuffer.resize(neighs.size());
		Communication::receiveUpperUnknownSize(sBuffer, rBuffer, neighs, group);

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
	pfile.offsets = Communication::getDistribution<size_t>(pfile.end - pfile.begin, group);
}

