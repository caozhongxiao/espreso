
#include "communication.h"

#include "../../config/ecf/input/loaderconfiguration.h"
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

MPITools::Operations::Operations()
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

MPITools::MPIGroup::MPIGroup()
{
	MPI_Comm_dup(environment->MPICommunicator, &communicator);
	rank = environment->MPIrank;
	size = environment->MPIsize;
}

MPITools::MPICommunicator::MPICommunicator()
{
	int color, length;
	std::vector<char> name(MPI_MAX_PROCESSOR_NAME);
	MPI_Get_processor_name(name.data(), &length);

	std::vector<int> rCounts(environment->MPIsize), rDispl(environment->MPIsize), colors(environment->MPIsize);
	std::vector<char> names;
	MPI_Gather(&length, 1, MPI_INT, rCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	for (size_t i = 1; i < rCounts.size(); i++) {
		rDispl[i] += rDispl[i - 1] + rCounts[i - 1];
	}
	names.resize(rDispl.back() + rCounts.back());

	MPI_Gatherv(name.data(), length * sizeof(char), MPI_BYTE, names.data(), rCounts.data(), rDispl.data(), MPI_BYTE, 0, MPI_COMM_WORLD);

	std::map<std::string, size_t> nodes;
	for (int i = 0; i < environment->MPIsize; i++) {
		std::string str(names.begin() + rDispl[i], names.begin() + rDispl[i] + rCounts[i]);
		auto it = nodes.find(str);
		if (it == nodes.end()) {
			size_t s = nodes.size();
			nodes[str] = s;
		}
		colors[i] = nodes[str];
	}

	MPI_Scatter(colors.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Comm_split(environment->MPICommunicator, color, environment->MPIrank, &within.communicator);
	MPI_Comm_rank(within.communicator, &within.rank);
	MPI_Comm_size(within.communicator, &within.size);

	MPI_Comm_split(MPI_COMM_WORLD, within.rank, environment->MPIrank, &across.communicator);
	MPI_Comm_rank(across.communicator, &across.rank);
	MPI_Comm_size(across.communicator, &across.size);
}

void Communication::createCommunicator(const LoaderConfiguration &configuration, MPITools::MPICommunicator &communicator)
{
	int color;
	MPITools::MPIGroup *group;

	switch (configuration.load_by) {
	case LoaderConfiguration::LoadBy::NODES:
		group = &MPITools::nodes().across;
		break;
	case LoaderConfiguration::LoadBy::PROCESSES:
		group = &MPITools::procs();
		break;
	}

	switch (configuration.load_pattern) {
	case LoaderConfiguration::LoadPattern::PREFIX:
		color = group->rank < configuration.pattern_value;
		MPI_Comm_split(environment->MPICommunicator, color, environment->MPIrank, &communicator.within.communicator);
		MPI_Comm_rank(communicator.within.communicator, &communicator.within.rank);
		MPI_Comm_size(communicator.within.communicator, &communicator.within.size);

		MPI_Comm_split(MPI_COMM_WORLD, communicator.within.rank, environment->MPIrank, &communicator.across.communicator);
		MPI_Comm_rank(communicator.across.communicator, &communicator.across.rank);
		MPI_Comm_size(communicator.across.communicator, &communicator.across.size);
		break;
	case LoaderConfiguration::LoadPattern::SUBSET:
		color = group->rank / configuration.pattern_value;
		MPI_Comm_split(environment->MPICommunicator, color, environment->MPIrank, &communicator.within.communicator);
		MPI_Comm_rank(communicator.within.communicator, &communicator.within.rank);
		MPI_Comm_size(communicator.within.communicator, &communicator.within.size);

		MPI_Comm_split(MPI_COMM_WORLD, communicator.within.rank, environment->MPIrank, &communicator.across.communicator);
		MPI_Comm_rank(communicator.across.communicator, &communicator.across.rank);
		MPI_Comm_size(communicator.across.communicator, &communicator.across.size);
		break;
	}
}

