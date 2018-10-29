
#include "communication.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/processesreduction.h"
#include "../../mesh/store/statisticsstore.h"

#include <string>
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace espreso;

int MPISubset::nodeRank = -1;

template<typename Ttype>
static void _scan(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) += *(static_cast<Ttype*>(in) + i);
	}
}

static void _mergeStatistics(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Statistics);
	for (int i = 0; i < size; i++) {
		(static_cast<Statistics*>(out) + i)->min = std::min((static_cast<Statistics*>(out) + i)->min, (static_cast<Statistics*>(in) + i)->min);
		(static_cast<Statistics*>(out) + i)->max = std::max((static_cast<Statistics*>(out) + i)->max, (static_cast<Statistics*>(in) + i)->max);
		(static_cast<Statistics*>(out) + i)->avg += (static_cast<Statistics*>(in) + i)->avg;
		(static_cast<Statistics*>(out) + i)->norm += (static_cast<Statistics*>(in) + i)->norm;
	}
}

template<typename Ttype>
static void _max(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) = std::max(*(static_cast<Ttype*>(in) + i), *(static_cast<Ttype*>(out) + i));
	}
}

template<typename Ttype>
static void _min(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) = std::min(*(static_cast<Ttype*>(in) + i), *(static_cast<Ttype*>(out) + i));
	}
}

template<typename Ttype>
static void _sum(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) += *(static_cast<Ttype*>(in) + i);
	}
}

MPIOperations::MPIOperations()
{
	MPI_Op_create(_mergeStatistics, 1, &mergeStatistics);

	MPI_Op_create(_scan<size_t>, 1, &SIZET.scan);
	MPI_Op_create(_max<size_t>, 1, &SIZET.max);
	MPI_Op_create(_min<size_t>, 1, &SIZET.min);
	MPI_Op_create(_sum<size_t>, 1, &SIZET.sum);

	MPI_Op_create(_scan<int>, 1, &INT.scan);
	MPI_Op_create(_max<int>, 1, &INT.max);
	MPI_Op_create(_min<int>, 1, &INT.min);
	MPI_Op_create(_sum<int>, 1, &INT.sum);

	MPI_Op_create(_scan<long>, 1, &LONG.scan);
	MPI_Op_create(_max<long>, 1, &LONG.max);
	MPI_Op_create(_min<long>, 1, &LONG.min);
	MPI_Op_create(_sum<long>, 1, &LONG.sum);
}

MPIGroup::MPIGroup()
{
	MPI_Comm_dup(environment->MPICommunicator, &communicator);
	rank = environment->MPIrank;
	size = environment->MPIsize;
}

void MPISubset::fillNodeColor()
{
	int length, maxlength;
	std::vector<char> name(MPI_MAX_PROCESSOR_NAME);
	MPI_Get_processor_name(name.data(), &length);

	MPI_Allreduce(&length, &maxlength, 1, MPI_INT, MPI_MAX, environment->MPICommunicator);

	std::vector<char> names(environment->MPIsize * maxlength);

	MPI_Gather(name.data(), maxlength * sizeof(char), MPI_CHAR, names.data(), maxlength * sizeof(char), MPI_BYTE, 0, environment->MPICommunicator);

	std::vector<int> permutation(environment->MPIsize), colors(environment->MPIsize);

	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) {
		return std::lexicographical_compare(
				names.begin() + maxlength * i, names.begin() + maxlength * (i + 1),
				names.begin() + maxlength * j, names.begin() + maxlength * (j + 1));
	});
	for (size_t i = 1; i < colors.size(); i++) {
		if (memcmp(names.data() + maxlength * i, names.data() + maxlength * (i - 1), maxlength) ){
			colors[i] = colors[i - 1] + 1;
		} else {
			colors[i] = colors[i - 1];
		}
	}

	MPI_Scatter(colors.data(), 1, MPI_INT, &nodeRank, 1, MPI_INT, 0, environment->MPICommunicator);
}

MPISubset::MPISubset(const ProcessesReduction &reduction, MPIGroup &origin)
: origin(origin)
{
	int rank, color;

	switch (reduction.granularity) {
	case ProcessesReduction::Granularity::NODES:
		fillNodeColor();
		rank = nodeRank;
		break;
	case ProcessesReduction::Granularity::PROCESSES:
		rank = origin.rank;
		break;
	}

	switch (reduction.pattern) {
	case ProcessesReduction::Pattern::PREFIX:
		color = rank % reduction.reduction_ratio;
		break;
	case ProcessesReduction::Pattern::SUBSET:
		color = rank / reduction.reduction_ratio;
		break;
	}

	MPI_Comm_split(origin.communicator, color, origin.rank, &within.communicator);
	MPI_Comm_rank(within.communicator, &within.rank);
	MPI_Comm_size(within.communicator, &within.size);

	MPI_Comm_split(origin.communicator, within.rank, origin.rank, &across.communicator);
	MPI_Comm_rank(across.communicator, &across.rank);
	MPI_Comm_size(across.communicator, &across.size);
}

MPISubset& MPITools::nodes()
{
	ProcessesReduction reduction;
	reduction.granularity = ProcessesReduction::Granularity::NODES;
	reduction.pattern = ProcessesReduction::Pattern::SUBSET;
	reduction.reduction_ratio = 1;

	static MPISubset instance(reduction, MPITools::procs());
	return instance;
}

