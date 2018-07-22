
#include "communication.h"

#include "../../mesh/store/statisticsstore.h"

namespace espreso {

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

}

