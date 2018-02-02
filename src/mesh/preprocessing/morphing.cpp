
#include "meshpreprocessing.h"

#include "../mesh.h"
#include "../store/nodestore.h"
#include "../store/boundaryregionstore.h"
#include "../store/elementsregionstore.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/matrices/sparseCSRMatrix.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/logging/logging.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/meshmorphing.h"

using namespace espreso;

void MeshPreprocessing::morphRBF(const std::string &name, const RBFTargetConfiguration &configuration, int dimension)
{
	size_t threads = environment->OMP_NUM_THREADS;

	start("apply morphing '" + name + "'");

	if (_mesh->nodes->originCoordinates == NULL) {
		_mesh->nodes->originCoordinates = new serializededata<eslocal, Point>(*_mesh->nodes->coordinates);
	}

	std::vector<Point> sPoints, rPoints;
	std::vector<double> sDisplacement, rDisplacement;

	size_t nSize = 0;
	for(auto it = configuration.morphers.begin(); it != configuration.morphers.end(); ++it) {
		const BoundaryRegionStore *region = _mesh->bregion(it->first);
		nSize += region->nodes->datatarray().size();
	}

	sPoints.reserve(nSize);
	sDisplacement.reserve(dimension * nSize);

	for(auto it = configuration.morphers.begin(); it != configuration.morphers.end(); ++it) {
		const BoundaryRegionStore *region = _mesh->bregion(it->first);
		const auto &nodes = region->nodes->datatarray();
		const auto &coordinates = _mesh->nodes->coordinates->datatarray();

		switch (it->second.transformation) {
		case MORPHING_TRANSFORMATION::FIXED: {
			for (size_t i = 0; i < region->nintervals.size(); i++) {
				if (region->nintervals[i].sourceProcess == environment->MPIrank) {
					for (auto n = nodes.begin() + region->nintervals[i].begin; n != nodes.begin() + region->nintervals[i].end; ++n) {
						sPoints.push_back(coordinates[*n]);
					}
					sDisplacement.insert(sDisplacement.end(), dimension * (region->nintervals[i].end - region->nintervals[i].begin), 0);
				}
			}
		} break;

		case MORPHING_TRANSFORMATION::TRANSLATION: {
			for (size_t i = 0; i < region->nintervals.size(); i++) {
				if (region->nintervals[i].sourceProcess == environment->MPIrank) {
					size_t prevsize = sPoints.size();
					for (auto n = nodes.begin() + region->nintervals[i].begin; n != nodes.begin() + region->nintervals[i].end; ++n) {
						sPoints.push_back(coordinates[*n]);
					}
					sDisplacement.resize(sDisplacement.size() + dimension * (region->nintervals[i].end - region->nintervals[i].begin));
					it->second.translation.x.evaluator->evaluate(region->nintervals[i].end - region->nintervals[i].begin, dimension, sPoints.data() + prevsize, NULL, 0, sDisplacement.data() + dimension * prevsize + 0);
					it->second.translation.y.evaluator->evaluate(region->nintervals[i].end - region->nintervals[i].begin, dimension, sPoints.data() + prevsize, NULL, 0, sDisplacement.data() + dimension * prevsize + 1);
					if (dimension == 3) {
						it->second.translation.z.evaluator->evaluate(region->nintervals[i].end - region->nintervals[i].begin, dimension, sPoints.data() + prevsize, NULL, 0, sDisplacement.data() + dimension * prevsize + 2);
					}
				}
			}
		} break;

		case MORPHING_TRANSFORMATION::OFFSET:
		case MORPHING_TRANSFORMATION::SCALING:
		case MORPHING_TRANSFORMATION::ROTATION:
			// TODO:
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement mesh morphing tranformation.";
		}
	}

	if (!Communication::gatherUnknownSize(sPoints, rPoints)) {
		ESINFO(ERROR) << "ESPRESO internal error: gather morphed points";
	}

	if (!Communication::gatherUnknownSize(sDisplacement, rDisplacement)) {
		ESINFO(ERROR) << "ESPRESO internal error: gather morphed displacement";
	}

	std::vector<double> W, Q;

	if (environment->MPIrank == 0) {

		eslocal rowsFromCoordinates = rDisplacement.size() / dimension;

		// TODO: optimize
		DenseMatrix rhs(rowsFromCoordinates + dimension + 1, dimension);
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			for (int d = 0; d < dimension; d++) {
				rhs(r, d) = rDisplacement[r * dimension + d];
			}
		}

		DenseMatrix M(rowsFromCoordinates + dimension + 1, rowsFromCoordinates + dimension + 1);

		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M(rowsFromCoordinates + 0, r) = M(r, rowsFromCoordinates + 0) = rPoints[r].x;
			M(rowsFromCoordinates + 1, r) = M(r, rowsFromCoordinates + 1) = rPoints[r].y;
			if (dimension == 3) {
				M(rowsFromCoordinates + 2, r) = M(r, rowsFromCoordinates + 2) = rPoints[r].z;
			}

			M(rowsFromCoordinates + dimension, r) = M(r, rowsFromCoordinates + dimension) = 1;

			for(int rr = 0; rr < r; rr++) {
				M(r, rr) = M(rr, r) = configuration.function.evaluator->evaluate((rPoints[r] - rPoints[rr]).length());
			}
		}

		switch (configuration.solver) {
		case MORPHING_RBF_SOLVER::ITERATIVE: {
			SparseCSRMatrix<eslocal> MSparse(M);
			rhs.transpose();

			DenseMatrix wq(rhs.rows(), rhs.columns());
			for(eslocal r = 0 ; r < rhs.rows(); r++) {
				MSparse.gmresSolve(&rhs(r, 0), &wq(r, 0), configuration.solver_precision, 600);
			}

			for(eslocal c = 0; c < wq.columns(); c++) {
				for(eslocal r = 0; r < wq.rows(); r++) {
					if (c < rowsFromCoordinates) {
						W.push_back(wq(r, c));
					}else {
						Q.push_back(wq(r, c));
					}
				}
			}
		} break;

		case MORPHING_RBF_SOLVER::DIRECT:
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement mesh morphing solver.";
		}
	}

	if (!Communication::broadcastUnknownSize(rPoints)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast points.";
	}
	if (!Communication::broadcastUnknownSize(W)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast W.";
	}
	if (!Communication::broadcastUnknownSize(Q)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast Q.";
	}

	ElementsRegionStore *tregion = _mesh->eregion(configuration.target);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = tregion->nodes->begin(t)->begin(); n != tregion->nodes->end(t)->begin(); ++n) {
			Point &morphed = _mesh->nodes->coordinates->datatarray()[*n];
			Point origin = morphed;

			for (size_t i = 0; i < rPoints.size(); i++) {
				double R = configuration.function.evaluator->evaluate((rPoints[i] - origin).length());

				morphed.x += R * W[i * dimension + 0];
				morphed.y += R * W[i * dimension + 1];
				if (dimension == 3) {
					morphed.z += R * W[i * dimension + 2];
				}
			}
			if (dimension == 3) {
				morphed.x += origin.x * Q[0 * dimension + 0];
				morphed.x += origin.y * Q[1 * dimension + 0];
				morphed.x += origin.z * Q[2 * dimension + 0];
				morphed.x +=            Q[3 * dimension + 0];

				morphed.y += origin.x * Q[0 * dimension + 1];
				morphed.y += origin.y * Q[1 * dimension + 1];
				morphed.y += origin.z * Q[2 * dimension + 1];
				morphed.y +=            Q[3 * dimension + 1];

				morphed.z += origin.x * Q[0 * dimension + 2];
				morphed.z += origin.y * Q[1 * dimension + 2];
				morphed.z += origin.z * Q[2 * dimension + 2];
				morphed.z +=            Q[3 * dimension + 2];
			}
			if (dimension == 2) {
				morphed.x += origin.x * Q[0 * dimension + 0];
				morphed.x += origin.y * Q[1 * dimension + 0];
				morphed.x +=            Q[2 * dimension + 0];

				morphed.y += origin.x * Q[0 * dimension + 1];
				morphed.y += origin.y * Q[1 * dimension + 1];
				morphed.y +=            Q[2 * dimension + 1];
			}
		}
	}

	if (_morphing == NULL) {
		_morphing = _mesh->nodes->appendData(3, { "RBF_MORPHING" }, true);
		_morphing->gatheredData.resize(3 * _mesh->nodes->uniqueSize);
	}

	#pragma omp parallel for
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); i++) {
		const auto &origin = _mesh->nodes->originCoordinates->datatarray();
		const auto &morphed = _mesh->nodes->coordinates->datatarray();
		if (_mesh->nodes->pintervals[i].sourceProcess == environment->MPIrank) {
			eslocal offset = _mesh->nodes->pintervals[i].globalOffset - _mesh->nodes->uniqueOffset;
			for (eslocal n = _mesh->nodes->pintervals[i].begin; n < _mesh->nodes->pintervals[i].end; ++n, ++offset) {
				_morphing->gatheredData[3 * offset + 0] = (morphed[n] - origin[n]).x;
				_morphing->gatheredData[3 * offset + 1] = (morphed[n] - origin[n]).y;
				_morphing->gatheredData[3 * offset + 2] = (morphed[n] - origin[n]).z;
			}
		}
	}

	finish("apply morphing '" + name + "'");
}

