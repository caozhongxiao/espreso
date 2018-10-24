
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

#include "../../config/reader/tokenizer.h"

#include "mpi.h"

using namespace espreso;

void MeshPreprocessing::processMorpher(const RBFTargetTransformationConfiguration &target, int dimension,
		std::vector<Point> &sPoints, eslocal startPoint, std::vector<double> &sDisplacement) {

	int pointsToProcess = sPoints.size() - startPoint;

	switch (target.transformation) {
		case MORPHING_TRANSFORMATION::FIXED: {
			sDisplacement.insert(sDisplacement.end(), dimension * sPoints.size(), 0);
		} break;

		case MORPHING_TRANSFORMATION::TRANSLATION: {

			sDisplacement.resize(sDisplacement.size() + dimension * pointsToProcess);
			target.translation.x.evaluator->evaluate(pointsToProcess, dimension, sPoints.data() + startPoint, NULL, 0, sDisplacement.data() + dimension * startPoint + 0);
			target.translation.y.evaluator->evaluate(pointsToProcess, dimension, sPoints.data() + startPoint, NULL, 0, sDisplacement.data() + dimension * startPoint + 1);
			if (dimension == 3) {
				target.translation.z.evaluator->evaluate(pointsToProcess, dimension, sPoints.data() + startPoint, NULL, 0, sDisplacement.data() + dimension * startPoint + 2);
			}
		} break;

		case MORPHING_TRANSFORMATION::OFFSET: {
		} break;

		case MORPHING_TRANSFORMATION::SCALING: {

			std::vector<double> result;
			target.coordinate_system.createTranslationMatrixToCenter(result);
			std::vector<double> scaling;
			std::vector<double> sv(3);
			target.scaling.x.evaluator->evaluate(1, NULL, NULL, 0, &sv[0]);
			target.scaling.y.evaluator->evaluate(1, NULL, NULL, 0, &sv[1]);
			target.scaling.z.evaluator->evaluate(1, NULL, NULL, 0, &sv[2]);
			target.coordinate_system.createScalingMatrix(scaling, sv[0]/100.0, sv[1]/100.0, sv[2]/100.0);
			std::vector<double> TtoZero;
			target.coordinate_system.createTranslationMatrixToZero(TtoZero);
			target.coordinate_system.multiplyTransformationMatrices(scaling, result);
			target.coordinate_system.multiplyTransformationMatrices(TtoZero, result);

			for(auto it = sPoints.begin()+startPoint; it!=sPoints.end(); it++) {
				const Point& p1 = *it;
				Point p2= target.coordinate_system.applyTransformation(result, p1);
				//std::cout<<"P1 "<<p1<<" P2 "<<p2<<"\n";
				sDisplacement.push_back(p2.x - p1.x);
				sDisplacement.push_back(p2.y - p1.y);
				if (dimension == 3) {
					sDisplacement.push_back(p2.z - p1.z);
				}
			}
		} break;

		case MORPHING_TRANSFORMATION::ROTATION: {

			std::vector<double> result;
			target.coordinate_system.createTranslationMatrixToCenter(result);
			std::vector<double> rotation;
			target.coordinate_system.createRotationMatrix(rotation);
			std::vector<double> TtoZero;
			target.coordinate_system.createTranslationMatrixToZero(TtoZero);
			target.coordinate_system.multiplyTransformationMatrices(rotation, result);
			target.coordinate_system.multiplyTransformationMatrices(TtoZero, result);

			for(auto it = sPoints.begin()+startPoint; it!=sPoints.end(); it++) {
				const Point& p1 = *it;
				Point p2= target.coordinate_system.applyTransformation(result, p1);
				//std::cout<<"P1 "<<p1<<" P2 "<<p2<<"\n";
				sDisplacement.push_back(p2.x - p1.x);
				sDisplacement.push_back(p2.y - p1.y);
				if (dimension == 3) {
					sDisplacement.push_back(p2.z - p1.z);
				}
			}
		}break;
	 	 default:
	 		 ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement mesh morphing tranformation.";
	}
}

eslocal MeshPreprocessing::prepareMatrixM(std::vector<Point> &rPoints,
		std::vector<double> &rDisplacement,
		int dimension, const RBFTargetConfiguration &configuration,
		std::vector<double> &M_values,
		bool use_x, bool use_y, bool use_z
) {

	eslocal rowsFromCoordinates = rPoints.size();
	eslocal realsize = rowsFromCoordinates;

	M_values.clear();

	for (eslocal r = 0; r < rowsFromCoordinates; r++) {
		for(eslocal rr = 0; rr <= r; rr++) {
			M_values.push_back(configuration.function.evaluator->evaluate((rPoints[r] - rPoints[rr]).length()));
		}
	}
	if (use_x) {
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].x);
		}
		M_values.push_back(0);
		realsize++;
	}

	if (use_y) {
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].y);
		}
		M_values.push_back(0);
		realsize++;
		if (use_x)	M_values.push_back(0);
	}

	if (dimension == 3 && use_z) {
		for (eslocal r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].z);
		}
		if (use_x) M_values.push_back(0);
		if (use_y) M_values.push_back(0);
		M_values.push_back(0);
		realsize++;
	}

	for (eslocal r = 0; r < rowsFromCoordinates; r++) {
		M_values.push_back(1);
	}
	if (use_x) M_values.push_back(0);
	if (use_y) M_values.push_back(0);
	if (dimension == 3 && use_z) M_values.push_back(0);
	M_values.push_back(0);
	realsize++;
	return realsize;
}

void MeshPreprocessing::readExternalFile(
		const RBFTargetConfiguration &configuration, int dimension,
		std::map<std::string, std::vector<Point>> &external_data) {

	Tokenizer tokenizer(configuration.external_ffd.path);
	Tokenizer::Token token;

	auto myNextToken = [] (Tokenizer &tokenizer) -> Tokenizer::Token {
		Tokenizer::Token token;
		do {
			token = tokenizer.next();
		}while(token==Tokenizer::Token::LINE_END);
		return token;
	};
	auto myExpectToken =
			[] (Tokenizer &tokenizer, Tokenizer::Token real, Tokenizer::Token expected, const std::string &expectedToken) -> void {
				if (real!=expected) {
					ESINFO(GLOBAL_ERROR) << "Line: " << tokenizer.line() << ", unexpected symbol, was expecting "<<expectedToken;
				}
			};
	auto myConvert =
			[] (Tokenizer &tokenizer) -> double {
				size_t size;
				double result = std::stod(tokenizer.value(), &size);
				if (size != tokenizer.value().size()) {
					ESINFO(GLOBAL_ERROR) << "Line: " << tokenizer.line() << ", error while converting "<<tokenizer.value()<<" into a number.";
				}
				return result;
			};

	token = myNextToken(tokenizer);
	while (token != Tokenizer::Token::END) {
		myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
				" a region name.");
		std::string name = tokenizer.value();
		if (external_data.find(name) != external_data.end()) {
			ESINFO(GLOBAL_ERROR) << "Line: " << tokenizer.line()
					<< ", region " + tokenizer.value()
					<< " defined multiple times in the external file.";
		}
		token = myNextToken(tokenizer);
		myExpectToken(tokenizer, token, Tokenizer::Token::OBJECT_OPEN,
				" symbol \"{\".");
		token = myNextToken(tokenizer);
		while (token != Tokenizer::Token::OBJECT_CLOSE) {
			Point p;
			myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
					" a number.");
			p.x = myConvert(tokenizer);
			token = myNextToken(tokenizer);
			myExpectToken(tokenizer, token, Tokenizer::Token::DELIMITER,
					" symbol \",\".");

			token = myNextToken(tokenizer);
			myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
					" a number.");
			p.y = myConvert(tokenizer);

			token = myNextToken(tokenizer);

			if (dimension == 3) {
				myExpectToken(tokenizer, token, Tokenizer::Token::DELIMITER,
						" symbol \",\".");

				token = myNextToken(tokenizer);
				myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
						" a number.");
				p.z = myConvert(tokenizer);

				token = myNextToken(tokenizer);
			}
			myExpectToken(tokenizer, token, Tokenizer::Token::EXPRESSION_END,
					" symbol \";\".");
			token = myNextToken(tokenizer);

			external_data[name].push_back(p);
		}
		token = myNextToken(tokenizer);
	}
}

void MeshPreprocessing::morphRBF(const std::string &name, const RBFTargetConfiguration &configuration, int dimension)
{
	size_t threads = environment->OMP_NUM_THREADS;
	MATH::setNumberOfThreads(threads);

	start("processing morphing '" + name + "'");

	start("preparing data for morphing '" + name + "'");

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

		size_t prevsize = sPoints.size();
		for (size_t i = 0; i < region->nintervals.size(); i++) {
			if (region->nintervals[i].sourceProcess == environment->MPIrank) {
				for (auto n = nodes.begin() + region->nintervals[i].begin; n != nodes.begin() + region->nintervals[i].end; ++n) {
					sPoints.push_back(coordinates[*n]);
				}
			}
		}
		processMorpher(it->second, dimension, sPoints, prevsize, sDisplacement);
	}

	if (environment->MPIrank == 0) {
		if (configuration.external_ffd.path!="") {
			std::map<std::string, std::vector<Point>> external_data;

			readExternalFile(configuration, dimension, external_data);

			for(auto it = configuration.external_ffd.morphers.begin(); it != configuration.external_ffd.morphers.end(); ++it) {
				size_t prevsize = sPoints.size();
				if (external_data.find(it->first)==external_data.end()) {
					ESINFO(GLOBAL_ERROR) << "Region " << it->first << " does not exist in external file: "<<configuration.external_ffd.path;
				}

				for (auto n = external_data[it->first].begin(); n != external_data[it->first].end(); ++n) {
					sPoints.push_back(*n);
				}
				processMorpher(it->second, dimension, sPoints, prevsize, sDisplacement);
			}
		}
	}
	finish("preparing data for morphing '" + name + "'");

	start("transmitting data for morphing '" + name + "'");
	if (!Communication::gatherUnknownSize(sPoints, rPoints)) {
		ESINFO(ERROR) << "ESPRESO internal error: gather morphed points";
	}

	if (!Communication::broadcastUnknownSize(rPoints)) {
			ESINFO(ERROR) << "ESPRESO internal error: broadcast points.";
	}

	if (!Communication::gatherUnknownSize(sDisplacement, rDisplacement)) {
		ESINFO(ERROR) << "ESPRESO internal error: gather morphed displacement";
	}

	finish("transmitting data for morphing '" + name + "'");

	start("solving data for morphing '" + name + "'");
	std::vector<double> wq_values;

	if (environment->MPIrank == 0 ||
			(configuration.solver == MORPHING_RBF_SOLVER::ITERATIVE &&
			 environment->MPIrank < dimension &&
			 environment->MPIsize >= dimension)) {

		std::vector<double> M_values;
		//eslocal rowsFromCoordinates = rPoints.size();
		size_t M_size = rPoints.size() + dimension + 1;

		size_t realSize = prepareMatrixM(rPoints, rDisplacement, dimension, configuration, M_values);

		if (realSize != M_size) {
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: error while building matrix M.";
		}

		switch (configuration.solver) {
		case MORPHING_RBF_SOLVER::ITERATIVE: {

			std::vector<double> rhs_values;

			if (environment->MPIrank == 0) {

				rhs_values.resize(M_size*dimension);
				wq_values.resize(M_size*dimension);

				for (int d = 0; d < dimension; d++) {
					size_t r;
					for (r = 0; r < rPoints.size(); r++) {
						rhs_values[M_size*d + r] = rDisplacement[r * dimension + d];
					}
					for ( ; r < M_size; r++) {
						rhs_values[M_size*d + r] = 0;
					}
				}
			}

			std::vector<eslocal> iterations;

			if (environment->MPIrank == 0 && environment->MPIsize < dimension) {
				for(eslocal d = 0 ; d < dimension; d++) {
					eslocal itercount;
					MATH::SOLVER::GMRESUpperSymetricColumnMajorMat(
							M_size, &M_values[0],
							&rhs_values[d * M_size], &wq_values[d * M_size],
							configuration.solver_precision, configuration.solver_max_iter, itercount);
					iterations.push_back(itercount);
				}
			}else {

				if (environment->MPIrank == 0) {

					for(eslocal d = 1 ; d < dimension; d++) {
						MPI_Send(&rhs_values[M_size*d], M_size, MPI_DOUBLE,
									d, 0, environment->MPICommunicator);
					}

				}else{

					rhs_values.resize(M_size);
					wq_values.resize(M_size);

					MPI_Recv(rhs_values.data(),M_size,
							MPI_DOUBLE, 0, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);

				}
				eslocal itercount;
				MATH::SOLVER::GMRESUpperSymetricColumnMajorMat(
						M_size, &M_values[0],
						&rhs_values[0], &wq_values[0],
						configuration.solver_precision, configuration.solver_max_iter, itercount);

				iterations.push_back(itercount);

				if (environment->MPIrank == 0) {

					for (int d = 1; d < dimension; d++) {
						MPI_Recv(&wq_values[d*M_size],M_size, MPI_DOUBLE,
								d, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
						MPI_Recv(&itercount,1, MPI_INT,
										d, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
						iterations.push_back(itercount);
					}

				}else {
					MPI_Send(wq_values.data(), M_size, MPI_DOUBLE,
							0, 0, environment->MPICommunicator);
					MPI_Send(&itercount, 1, MPI_INT,
								0, 0, environment->MPICommunicator);
				}

			}

		} break;

		case MORPHING_RBF_SOLVER::DIRECT: {

			wq_values.resize(M_size*dimension);

			for (int d = 0; d < dimension; d++) {
				size_t r;
				for (r = 0; r < rPoints.size(); r++) {
					wq_values[M_size*d + r] = rDisplacement[r * dimension + d];
				}
				for ( ; r < M_size; r++) {
					wq_values[M_size*d + r] = 0;
				}
			}

			/*std::cout<<"Points: "<<rPoints.size()<<" - "<<rPoints<<"\n";

			int count=0;
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

				bool use_x = true,use_y = true,use_z = true;

				std::vector<double> xs(rPoints.size()), ys(rPoints.size()), zs(rPoints.size());
				for(size_t i=0; i<rPoints.size(); i++) {
					Point &p = rPoints[i];
					xs[i] = p.x;
					ys[i] = p.y;
					if (dimension == 3) {
						zs[i] = p.z;
					}
				}

				auto checkAllSame = [] (std::vector<double> &points) -> bool {
					if (points.size()==0) return false;
					for(auto it = points.begin()+1; it!=points.end();++it) {
						if (points[0]!=*it) {
							return false;
						}
					}
					return true;
				};

				if (checkAllSame(xs)) {
					use_x=false;
				}
				if (checkAllSame(ys)) {
					use_y=false;
				}
				if (dimension == 3) {
					if (checkAllSame(zs)) {
						use_z=false;
					}
				}

				eslocal realSize = prepareMatrixM(rPoints, rDisplacement, dimension, configuration, M_values,
						use_x, use_y, use_z);

				wq_values.clear();
				wq_values.resize(realSize*dimension);

				for (int d = 0; d < dimension; d++) {
					size_t r;
					for (r = 0; r < rPoints.size(); r++) {
						wq_values[realSize*d + r] = rDisplacement[r * dimension + d];
					}
					for ( ; r < M_size; r++) {
						wq_values[realSize*d + r] = 0;
					}
				}

				/*std::cout << "Points: " << rPoints.size() << " - " << rPoints
						<< "\n";

				int count = 0;
				for (int i = 0; i < realSize; i++) {
					for (int j = 0; j <= i; j++) {
						printf("%f ", M_values[count]);
						count++;
					}
					printf("\n");
				}

				printf("WQ\n");
				for(int d=0;d<dimension;d++) {
					for(int i=0;i<wq_values.size()/dimension;i++) {
						printf("%f ", wq_values[wq_values.size()/dimension*d+i]);
					}
					printf("\n");
				}*/

				int result= MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
						realSize,  &M_values[0],
						dimension, &wq_values[0]);

				if (result!=0) {
					ESINFO(ERROR) << "ESPRESO error: Dense solver is unable to solve the system. Try iterative solver.";
				}

				auto insertRowInWQ =
						[&wq_values,dimension] (eslocal row) -> void {
							eslocal size = wq_values.size()/dimension;
							wq_values.insert(wq_values.begin()+row,0);
							wq_values.insert(wq_values.begin()+ size +1 + row,0);
							if (dimension ==3) {
								wq_values.insert(wq_values.begin()+ (size +1)*2 + row,0);
							}
						};

				if (!use_x) insertRowInWQ(rPoints.size());
				if (!use_y) insertRowInWQ(rPoints.size()+1);
				if (dimension == 3 && !use_z) insertRowInWQ(rPoints.size()+2);

			}
		} break;

		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement mesh morphing solver.";
		}
	}

	finish("solving data for morphing '" + name + "'");

	start("transmitting results for morphing '" + name + "'");
	if (!Communication::broadcastUnknownSize(wq_values)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast WQ.";
	}
	finish("transmitting results for morphing '" + name + "'");


	start("applying morphing '" + name + "'");

	eslocal wq_points = wq_values.size()/dimension;
	eslocal points_size = rPoints.size();

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

		}
	}
	finish("applying morphing '" + name + "'");


	if (_morphing == NULL) {
		_morphing = _mesh->nodes->appendData(3, { "RBF_MORPHING", "RBF_MORPHING_X", "RBF_MORPHING_Y", "RBF_MORPHING_Z" }, true);
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

	finish("processing morphing '" + name + "'");

	MATH::setNumberOfThreads(1);
}

