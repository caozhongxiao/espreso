
#include "communication.h"

#include "../../mesh/store/statisticsstore.h"

namespace espreso {

template<typename Ttype>
static void _sizeToOffsets(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	*(static_cast<Ttype*>(out)) += *(static_cast<Ttype*>(in));
}

static void _mergeStatistics(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	for (int i = 0; i < *len; i++) {
		static_cast<Statistics*>(out)->min = std::min(static_cast<Statistics*>(out)->min, static_cast<Statistics*>(in)->min);
		static_cast<Statistics*>(out)->max = std::max(static_cast<Statistics*>(out)->max, static_cast<Statistics*>(in)->max);
		static_cast<Statistics*>(out)->avg += static_cast<Statistics*>(in)->avg;
		static_cast<Statistics*>(out)->norm += static_cast<Statistics*>(in)->norm;
	}
}

MPITools::Operations::Operations()
{
	MPI_Op_create(_sizeToOffsets<eslocal>, 1, &sizeToOffsets);
	MPI_Op_create(_mergeStatistics, 1, &mergeStatistics);
}

}

