
#include "zones.h"

#include "input/openfoam/openfoam.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/communication.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

using namespace espreso;

OpenFOAMZones::OpenFOAMZones(ParallelFile &pfile)
: OpenFOAMCollectiveParser(pfile.begin, pfile.end), _pfile(pfile)
{

}

int OpenFOAMZones::getZones()
{
	int zones = 0, root = 0, rank = 0;
	const char *c = begin;
	if (*(c - 1) == '(') {
		c -= 3;
		while (*c != '\n') { c--; } // go before number of boundaries
		rank = info::mpi::MPIrank;
		zones = readInteger(c);
	}

	MPI_Allreduce(&rank, &root, 1, MPI_INT, MPI_SUM, info::mpi::MPICommunicator);
	MPI_Bcast(&zones, 1, MPI_INT, root, info::mpi::MPICommunicator);
	return zones;
}

void OpenFOAMZones::synchronize(int zones, std::vector<char> &names, std::vector<size_t> &offsets)
{
	names.resize(zones * 80);

	std::vector<char> mynames(zones * 80);
	memset(mynames.data(), '\0', mynames.size());
	int mybrackets = 0, scannedBrackets = 0;
	int myoffsets = 0, scannedOffsets = 0;
	std::vector<const char*> brackets;

	if (begin != end) {
		const char *c = begin + 1, *line = begin;
		while (c != end) {
			if (*c == '\n') {
				line = c + 1;
			}
			if (*c == '{') {
				if (!StringCompare::caseInsensitiveEq(readString(line), "flipMap")) {
					brackets.push_back(c);
				}
			}
			if (*c == '(' || *c == ')') {
				offsets.push_back(_pfile.offsets[info::mpi::MPIrank] + (c - _pfile.begin));
			}
			++c;
		}
		if (*c == '{') { // if '{' is the first of the next process, I have to parse the name
			brackets.push_back(c);
		}
	}

	mybrackets = brackets.size();
	MPI_Exscan(&mybrackets, &scannedBrackets, 1, MPI_INT, MPI_SUM, info::mpi::MPICommunicator);
	myoffsets = offsets.size();
	MPI_Exscan(&myoffsets, &scannedOffsets, 1, MPI_INT, MPI_SUM, info::mpi::MPICommunicator);
	std::vector<size_t> _offset(scannedOffsets);
	_offset.insert(_offset.end(), offsets.begin(), offsets.end());
	_offset.resize(2 * zones);
	offsets.resize(2 * zones);
	MPI_Allreduce(_offset.data(), offsets.data(), _offset.size() * sizeof(size_t), MPI_BYTE, MPITools::sizetOperations().sum, info::mpi::MPICommunicator);

	for (size_t i = 0; i < brackets.size(); i++) {
		const char *name = brackets[i];
		while (*name-- != '\n');
		while (*name != '\n') { name--; }
		++name;
		std::string zonename = readString(name);
		memcpy(mynames.data() + 80 * (scannedBrackets + i), zonename.data(), zonename.size());
	}

	MPI_Allreduce(mynames.data(), names.data(), names.size(), MPI_CHAR, MPI_SUM, info::mpi::MPICommunicator);
}

void OpenFOAMZones::readData(std::vector<esint> &indices, size_t begin, size_t end)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	if (begin + 1 > _pfile.offsets[info::mpi::MPIrank + 1]) {
		return;
	}
	if (end - 1 < _pfile.offsets[info::mpi::MPIrank]) {
		return;
	}

	begin = std::max(begin + 1, _pfile.offsets[info::mpi::MPIrank]);
	end = std::min(end - 1, _pfile.offsets[info::mpi::MPIrank + 1]);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);

	std::vector<std::vector<esint> > data(threads);
	size_t offset = begin - _pfile.offsets[info::mpi::MPIrank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;

		const char *c = _pfile.begin + offset + tdistribution[t];
		if (c > _pfile.begin) {
			while (c < _pfile.end && *(c - 1) != '\n') { ++c; }
		}
		while (c < _pfile.begin + offset + tdistribution[t + 1]) {
			tdata.push_back(readInteger(c));
			c += 1; // skip '\n'
		}

		data[t].swap(tdata);
	}

	for (size_t t = 0; t < threads; t++) {
		indices.insert(indices.end(), data[t].begin(), data[t].end());
	}
	std::sort(indices.begin(), indices.end());
}

bool OpenFOAMZones::readPoints(PlainOpenFOAMData &data)
{
	int zones = getZones();
	if (zones == 0) {
		return true;
	}

	std::vector<char> names;
	std::vector<size_t> offsets;
	synchronize(zones, names, offsets);

	for (int i = 0; i < zones; i++) {
		std::string name(names.data() + 80 * i);
		auto &indices = data.nregions[name];
		readData(indices, offsets[2 * i], offsets[2 * i + 1]);
	}

	return true;
}

bool OpenFOAMZones::readFaces(PlainOpenFOAMData &data)
{
	int zones = getZones();
	if (zones == 0) {
		return true;
	}

	std::vector<char> names;
	std::vector<size_t> offsets;
	synchronize(zones, names, offsets);

	for (int i = 0; i < zones; i++) {
		std::string name(names.data() + 80 * i);
		auto &indices = data.eregions[data.boundaryprefix + name];
		readData(indices, offsets[2 * i], offsets[2 * i + 1]);
	}

	return true;
}

bool OpenFOAMZones::readCells(PlainOpenFOAMData &data)
{
	int zones = getZones();
	if (zones == 0) {
		return true;
	}

	std::vector<char> names;
	std::vector<size_t> offsets;
	synchronize(zones, names, offsets);

	for (int i = 0; i < zones; i++) {
		std::string name(names.data() + 80 * i);
		auto &indices = data.eregions[data.elementprefix + name];
		readData(indices, offsets[2 * i], offsets[2 * i + 1]);
	}

	return true;
}






