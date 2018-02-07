
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

#include "../../wrappers/math/math.h"

using namespace espreso;

void MeshPreprocessing::morphRBF(const std::string &name, const RBFTargetConfiguration &configuration, int dimension)
{
	size_t threads = environment->OMP_NUM_THREADS;

	start("apply morphing '" + name + "'");

	ESINFO(OVERVIEW)<<"Processing morphing: "<<name;

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

			case MORPHING_TRANSFORMATION::OFFSET: {

			} break;

			case MORPHING_TRANSFORMATION::SCALING: {

				std::vector<double> result;
				it->second.coordinate_system.createTranslationMatrixToCenter(result);
				std::vector<double> scaling;
				std::vector<double> sv(3);
				it->second.scaling.x.evaluator->evaluate(1, NULL, NULL, 0, &sv[0]);
				it->second.scaling.y.evaluator->evaluate(1, NULL, NULL, 0, &sv[1]);
				it->second.scaling.z.evaluator->evaluate(1, NULL, NULL, 0, &sv[2]);
				it->second.coordinate_system.createScalingMatrix(scaling, sv[0]/100.0, sv[1]/100.0, sv[2]/100.0);
				std::vector<double> TtoZero;
				it->second.coordinate_system.createTranslationMatrixToZero(TtoZero);
				it->second.coordinate_system.multiplyTransformationMatrices(scaling, result);
				it->second.coordinate_system.multiplyTransformationMatrices(TtoZero, result);


				for (size_t i = 0; i < region->nintervals.size(); i++) {
					if (region->nintervals[i].sourceProcess == environment->MPIrank) {
						size_t prevsize = sPoints.size();
						for (auto n = nodes.begin() + region->nintervals[i].begin; n != nodes.begin() + region->nintervals[i].end; ++n) {
							const Point& p1 = coordinates[*n];
							sPoints.push_back(p1);
							Point p2= it->second.coordinate_system.applyTransformation(result, p1);
							//std::cout<<"P1 "<<p1<<" P2 "<<p2<<"\n";
							sDisplacement.push_back(p2.x - p1.x);
							sDisplacement.push_back(p2.y - p1.y);
							if (dimension == 3) {
								sDisplacement.push_back(p2.z - p1.z);
							}
						}
					}
				}
			} break;

			case MORPHING_TRANSFORMATION::ROTATION: {

				std::vector<double> result;
				it->second.coordinate_system.createTranslationMatrixToCenter(result);
				std::vector<double> rotation;
				it->second.coordinate_system.createRotationMatrix(rotation);
				std::vector<double> TtoZero;
				it->second.coordinate_system.createTranslationMatrixToZero(TtoZero);
				it->second.coordinate_system.multiplyTransformationMatrices(rotation, result);
				it->second.coordinate_system.multiplyTransformationMatrices(TtoZero, result);

				for (size_t i = 0; i < region->nintervals.size(); i++) {
					if (region->nintervals[i].sourceProcess == environment->MPIrank) {
						size_t prevsize = sPoints.size();
						for (auto n = nodes.begin() + region->nintervals[i].begin; n != nodes.begin() + region->nintervals[i].end; ++n) {
							const Point& p1 = coordinates[*n];
							sPoints.push_back(p1);
							Point p2= it->second.coordinate_system.applyTransformation(result, p1);
							//std::cout<<"P1 "<<p1<<" P2 "<<p2<<"\n";
							sDisplacement.push_back(p2.x - p1.x);
							sDisplacement.push_back(p2.y - p1.y);
							if (dimension == 3) {
								sDisplacement.push_back(p2.z - p1.z);
							}
						}
					}
				}

			}break;

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

	std::vector<double> wq_values;

	if (environment->MPIrank == 0) {

		eslocal rowsFromCoordinates = rDisplacement.size() / dimension;
		int M_size = rowsFromCoordinates + dimension + 1;

		/*DenseMatrix M(rowsFromCoordinates + dimension + 1, rowsFromCoordinates + dimension + 1);

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
		std::cout<<M;*/

		/*std::vector<double> M_values;//(M_size*(M_size+1)/2);
		for(eslocal i = 0;i<M_size;i++) {
			for(eslocal j = 0;j<=i;j++) {
				M_values.push_back(M(i,j));
			}
		}*/

		std::vector<double> M_values;

		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			for(eslocal rr = 0; rr <= r; rr++) {
				M_values.push_back(configuration.function.evaluator->evaluate((rPoints[r] - rPoints[rr]).length()));
			}
		}
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].x);
		}
		M_values.push_back(0);
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].y);
		}
		M_values.push_back(0);
		M_values.push_back(0);
		if (dimension == 3) {
			for (eslocal r = 0; r < rowsFromCoordinates; r++) {
				M_values.push_back(rPoints[r].z);
			}
			M_values.push_back(0);
			M_values.push_back(0);
			M_values.push_back(0);
		}
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(1);
		}
		M_values.push_back(0);
		M_values.push_back(0);
		M_values.push_back(0);
		M_values.push_back(0);


		switch (configuration.solver) {
		case MORPHING_RBF_SOLVER::ITERATIVE: {

			wq_values.resize(M_size*dimension);

			std::vector<double> rhs_values(M_size*dimension);
			for (int d = 0; d < dimension; d++) {
				eslocal r;
				for (r = 0; r < rowsFromCoordinates; r++) {
					rhs_values[M_size*d + r] = rDisplacement[r * dimension + d];
				}
				for ( ; r < M_size; r++) {
					rhs_values[M_size*d + r] = 0;
				}
			}

			for(eslocal d = 0 ; d < dimension; d++) {
				MATH::SOLVER::GMRESUpperSymetricColumnMajorMat(
					 M_size, &M_values[0],
					&rhs_values[d * M_size], &wq_values[d * M_size],
					configuration.solver_precision, 600);
			}

		} break;

		case MORPHING_RBF_SOLVER::DENSE: {

			wq_values.resize(M_size*dimension);

			for (int d = 0; d < dimension; d++) {
				eslocal r;
				for (r = 0; r < rowsFromCoordinates; r++) {
					wq_values[M_size*d + r] = rDisplacement[r * dimension + d];
				}
				for ( ; r < M_size; r++) {
					wq_values[M_size*d + r] = 0;
				}
			}

			/*int count=0;
			for(int i=0;i<M_size;i++) {
				for(int j=0;j<=i;j++) {
					printf("%f ", M_values[count]);
					count++;
				}
				printf("\n");
			}*/

			int result= MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
					M_size,  &M_values[0],
					dimension, &wq_values[0]);


			if (result!=0) {
				ESINFO(ERROR) << "ESPRESO error: Dense solver is unable to solve the system. Try iterative solver.";
			}
		} break;

		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement mesh morphing solver.";
		}
	}

	if (!Communication::broadcastUnknownSize(rPoints)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast points.";
	}
	if (!Communication::broadcastUnknownSize(wq_values)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast WQ.";
	}

	//std::vector<double> W, Q;
	eslocal wq_points = wq_values.size()/dimension;
	eslocal points_size = rPoints.size();

	/*for(eslocal c = 0; c < wq_points; c++) {
		for(eslocal r = 0; r < dimension; r++) {
			if (c < rPoints.size()) {
				W.push_back(wq_values[r * wq_points+ c]);
			}else {
				Q.push_back(wq_values[r * wq_points+ c]);
			}
		}
	}*/

	ElementsRegionStore *tregion = _mesh->eregion(configuration.target);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = tregion->nodes->begin(t)->begin(); n != tregion->nodes->end(t)->begin(); ++n) {
			Point &morphed = _mesh->nodes->coordinates->datatarray()[*n];
			Point origin = morphed;

			for (size_t i = 0; i < rPoints.size(); i++) {
				double R = configuration.function.evaluator->evaluate((rPoints[i] - origin).length());

				morphed.x += R * wq_values[i ];
				morphed.y += R * wq_values[i + wq_points];
				if (dimension == 3) {
					morphed.z += R * wq_values[i + wq_points * 2];
				}
			}
			if (dimension == 3) {
				morphed.x += origin.x * wq_values[points_size + 0];
				morphed.x += origin.y * wq_values[points_size + 1];
				morphed.x += origin.z * wq_values[points_size + 2];
				morphed.x +=            wq_values[points_size + 3];

				morphed.y += origin.x * wq_values[wq_points + points_size + 0];
				morphed.y += origin.y * wq_values[wq_points + points_size + 1];
				morphed.y += origin.z * wq_values[wq_points + points_size + 2];
				morphed.y +=            wq_values[wq_points + points_size + 3];

				morphed.z += origin.x * wq_values[wq_points*2 + points_size + 0];
				morphed.z += origin.y * wq_values[wq_points*2 + points_size + 1];
				morphed.z += origin.z * wq_values[wq_points*2 + points_size + 2];
				morphed.z +=            wq_values[wq_points*2 + points_size + 3];
			}

			if (dimension == 2) {
				morphed.x += origin.x * wq_values[points_size + 0];
				morphed.x += origin.y * wq_values[points_size + 1];
				morphed.x +=            wq_values[points_size + 2];

				morphed.y += origin.x * wq_values[wq_points + points_size + 0];
				morphed.y += origin.y * wq_values[wq_points + points_size + 1];
				morphed.y +=            wq_values[wq_points + points_size + 2];
			}

			/*for (size_t i = 0; i < rPoints.size(); i++) {
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
			}*/
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

