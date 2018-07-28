
#include "wparmetis.h"
#include "parmetis.h"

#include "../../basis/logging/logging.h"
#include "../../config/ecf/environment.h"

using namespace espreso;

eslocal ParMETIS::call(
			ParMETIS::METHOD method,
			eslocal *edistribution,
			eslocal *eframes, eslocal *eneighbors,
			eslocal dimensions, double *coordinates,
			eslocal verticesWeightCount, eslocal *verticesWeights, eslocal *edgeWeights,
			eslocal *partition)
{
	verticesWeightCount = std::max((eslocal)1, verticesWeightCount);

	eslocal wgtflag = 0;
	eslocal numflag = 0;
	eslocal parts = environment->MPIsize;
	std::vector<double> partFraction(verticesWeightCount * parts, 1.0 / parts);
	std::vector<double> unbalanceTolerance(verticesWeightCount, 1.02);
	eslocal options[4] = { 0, 0, 0, PARMETIS_PSR_UNCOUPLED };
	double itr = 1e6;
	eslocal edgecut;
	MPI_Comm communication = environment->MPICommunicator;

	if (verticesWeights != NULL) {
		wgtflag += 2;
	}
	if (edgeWeights != NULL) {
		wgtflag += 1;
	}

	switch (method) {

	case ParMETIS::METHOD::ParMETIS_V3_PartKway:
		if (coordinates != NULL) {
			if (METIS_OK != ParMETIS_V3_PartGeomKway(
					edistribution,
					eframes, eneighbors,
					verticesWeights, edgeWeights,
					&wgtflag, &numflag, &dimensions, coordinates, &verticesWeightCount,
					&parts, partFraction.data(), unbalanceTolerance.data(),
					options,
					&edgecut, partition,
					&communication)) {

				ESINFO(ERROR) << "PARMETIS_ERROR while partitiate mesh to MPI processes by KWay utilizing coordinates.";
			}
		} else {
			if (METIS_OK != ParMETIS_V3_PartKway(
					edistribution,
					eframes, eneighbors,
					verticesWeights, edgeWeights,
					&wgtflag, &numflag, &verticesWeightCount,
					&parts, partFraction.data(), unbalanceTolerance.data(),
					options,
					&edgecut, partition,
					&communication)) {

				ESINFO(ERROR) << "PARMETIS_ERROR while partitiate mesh to MPI processes by KWay utilizing coordinates.";
			}
		}
		break;

	case ParMETIS::METHOD::ParMETIS_V3_RefineKway:
		if (METIS_OK != ParMETIS_V3_RefineKway(
				edistribution,
				eframes, eneighbors,
				verticesWeights, edgeWeights,
				&wgtflag, &numflag, &verticesWeightCount,
				&parts, partFraction.data(), unbalanceTolerance.data(),
				options,
				&edgecut, partition,
				&communication)) {

			ESINFO(ERROR) << "PARMETIS_ERROR while refine mesh partition to MPI processes by KWay.";
		}
		break;

	case ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart:
		if (METIS_OK != ParMETIS_V3_AdaptiveRepart(
				edistribution,
				eframes, eneighbors,
				verticesWeights, NULL, edgeWeights,
				&wgtflag, &numflag, &verticesWeightCount,
				&parts, partFraction.data(), unbalanceTolerance.data(), &itr,
				options,
				&edgecut, partition,
				&communication)) {

			ESINFO(ERROR) << "PARMETIS_ERROR while adaptive repartition mesh to MPI processes by KWay.";
		}
		break;

	case ParMETIS::METHOD::ParMETIS_V3_PartGeom:
		if (coordinates == NULL) {
			ESINFO(ERROR) << "PARMETIS_ERROR:: PartGeom needs coordinates.";
		}
		if (METIS_OK != ParMETIS_V3_PartGeom(
				edistribution,
				&dimensions, coordinates,
				partition,
				&communication)) {

			ESINFO(ERROR) << "PARMETIS_ERROR while refine mesh partition to MPI processes by KWay.";
		}

	}

	return edgecut;
}


