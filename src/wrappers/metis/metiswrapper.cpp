
#include "metis.h"

#include "metiswrapper.h"
#include "basis/logging/logging.h"
#include "config/ecf/decomposition.h"

using namespace espreso;

esint METIS::call(
		const METISConfiguration &options,
		esint verticesCount,
		esint *eframes, esint *eneighbors,
		esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
		esint parts, esint *partition)
{

	verticesWeightCount = std::max((esint)1, verticesWeightCount);

	esint moptions[METIS_NOPTIONS];
	METIS_SetDefaultOptions(moptions);

	// HEURISTICS
	switch (options.objective_type) {
	case METISConfiguration::OBJECTIVE_TYPE::VOLUME:
		moptions[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_VOL;
		break;
	case METISConfiguration::OBJECTIVE_TYPE::EDGECUT:
		moptions[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
		break;
	default:
		ESINFO(ERROR) << "METIS WRAPPER error: invalid objective type.";
	}

	moptions[METIS_OPTION_CTYPE]     = METIS_CTYPE_SHEM;
	moptions[METIS_OPTION_IPTYPE]    = METIS_IPTYPE_EDGE;
	moptions[METIS_OPTION_RTYPE]     = METIS_RTYPE_FM;
	moptions[METIS_OPTION_NITER]     = 20;
	moptions[METIS_OPTION_UFACTOR]   = 100; // imbalance (1 + x) / 1000

	// MANDATORY
	moptions[METIS_OPTION_CONTIG]    = 1;
	moptions[METIS_OPTION_MINCONN]   = 0;
	moptions[METIS_OPTION_NUMBERING] = 0;
	moptions[METIS_OPTION_DBGLVL]    = 0;

	// KEPT DEFAULT
	// options[METIS_OPTION_NO2HOP]    = 0;
	// options[METIS_OPTION_NCUTS]     = 1;
	// options[METIS_OPTION_SEED]      = 0;

	esint edgecut = 0;

	if (parts > 1) {
		if (METIS_OK != METIS_PartGraphKway(
				&verticesCount, &verticesWeightCount,
				eframes, eneighbors,
				verticesWeights, NULL, edgeWeights,
				&parts, NULL, NULL,
				moptions, &edgecut, partition)) {

			ESINFO(ERROR) << "METIS_ERROR while KWay decomposition.";
		}
	} else {
		std::fill(partition, partition + verticesCount, 0);
	}

	return edgecut;
}


