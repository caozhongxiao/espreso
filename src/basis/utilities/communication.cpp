
#include "communication.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/processesreduction.h"
#include "../../mesh/store/statisticsstore.h"

#include <string>
#include <map>


using namespace espreso;

template<typename Ttype>
static void _scan(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	*(static_cast<Ttype*>(out)) += *(static_cast<Ttype*>(in));
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

MPISubset::MPISubset()
: origin(MPITools::procs())
{
/*	int color, length;
	std::vector<char> name(MPI_MAX_PROCESSOR_NAME);
	MPI_Get_processor_name(name.data(), &length);

	std::vector<int> rCounts(origin.size), rDispl(origin.size), colors(origin.size);
	std::vector<char> names;
	MPI_Gather(&length, 1, MPI_INT, rCounts.data(), 1, MPI_INT, 0, origin.communicator);

	for (size_t i = 1; i < rCounts.size(); i++) {
		rDispl[i] += rDispl[i - 1] + rCounts[i - 1];
	}
	names.resize(rDispl.back() + rCounts.back());

	MPI_Gatherv(name.data(), length * sizeof(char), MPI_BYTE, names.data(), rCounts.data(), rDispl.data(), MPI_BYTE, 0, origin.communicator);

	std::map<std::string, size_t> nodes;
	for (int i = 0; i < origin.size; i++) {
		std::string str(names.begin() + rDispl[i], names.begin() + rDispl[i] + rCounts[i]);
		auto it = nodes.find(str);
		if (it == nodes.end()) {
			size_t s = nodes.size();
			nodes[str] = s;
		}
		colors[i] = nodes[str];
	}

	MPI_Scatter(colors.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, origin.communicator);

	MPI_Comm_split(origin.communicator, color, origin.rank, &within.communicator);
	MPI_Comm_rank(within.communicator, &within.rank);
	MPI_Comm_size(within.communicator, &within.size);

	MPI_Comm_split(origin.communicator, within.rank, origin.rank, &across.communicator);
	MPI_Comm_rank(across.communicator, &across.rank);
	MPI_Comm_size(across.communicator, &across.size);*/
}

void Communication::createSubset(const ProcessesReduction &reduction, MPISubset &subset)
{
	int color;
//	MPIGroup *group;

//	switch (reduction.granularity) {
//	case ProcessesReduction::Granularity::NODES:
//		group = &MPITools::nodes().across;
//		break;
//	case ProcessesReduction::Granularity::PROCESSES:
//		group = &MPITools::procs();
//		break;
//	}

//	switch (reduction.pattern) {
//	case ProcessesReduction::Pattern::PREFIX:
//		color = group->rank % reduction.reduction_ratio;
//		break;
//	case ProcessesReduction::Pattern::SUBSET:
		color = MPITools::procs().rank / reduction.reduction_ratio;
//		break;
//	}

	MPI_Comm_split(subset.origin.communicator, color, subset.origin.rank, &subset.within.communicator);
	MPI_Comm_rank(subset.within.communicator, &subset.within.rank);
	MPI_Comm_size(subset.within.communicator, &subset.within.size);

	MPI_Comm_split(subset.origin.communicator, subset.within.rank, subset.origin.rank, &subset.across.communicator);
	MPI_Comm_rank(subset.across.communicator, &subset.across.rank);
	MPI_Comm_size(subset.across.communicator, &subset.across.size);
}

