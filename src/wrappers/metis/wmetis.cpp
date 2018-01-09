
#include "wmetis.h"
#include "metis.h"

#include "../../basis/logging/logging.h"
#include "../../config/ecf/environment.h"

#include <fstream>

using namespace espreso;

eslocal METIS::call(
		eslocal verticesCount,
		eslocal *eframes, eslocal *eneighbors,
		eslocal verticesWeightCount, eslocal *verticesWeights, eslocal *edgeWeights,
		eslocal parts, eslocal *partition)
{

	verticesWeightCount = std::max(1, verticesWeightCount);

	eslocal options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);

	// HEURISTICS
	options[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
	options[METIS_OPTION_CTYPE]     = METIS_CTYPE_SHEM;
	options[METIS_OPTION_IPTYPE]    = METIS_IPTYPE_EDGE;
	options[METIS_OPTION_RTYPE]     = METIS_RTYPE_FM;
	options[METIS_OPTION_NITER]     = 20;
	options[METIS_OPTION_UFACTOR]   = 100; // imbalance (1 + x) / 1000

	// MANDATORY
	options[METIS_OPTION_CONTIG]    = 1;
	options[METIS_OPTION_MINCONN]   = 0;
	options[METIS_OPTION_NUMBERING] = 0;
	options[METIS_OPTION_DBGLVL]    = 0;

	// KEPT DEFAULT
	// options[METIS_OPTION_NO2HOP]    = 0;
	// options[METIS_OPTION_NCUTS]     = 1;
	// options[METIS_OPTION_SEED]      = 0;

	eslocal edgecut = 0;

	if (parts > 1) {
		if (METIS_OK != METIS_PartGraphKway(
				&verticesCount, &verticesWeightCount,
				eframes, eneighbors,
				verticesWeights, NULL, edgeWeights,
				&parts, NULL, NULL,
				options, &edgecut, partition)) {

			ESINFO(ERROR) << "METIS_ERROR while KWay decomposition.";
		}
	} else {
		std::fill(partition, partition + verticesCount, 0);
	}

	return edgecut;
}


