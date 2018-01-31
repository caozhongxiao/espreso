
#include "meshpreprocessing.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/logging/logging.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/meshmorphing.h"

using namespace espreso;

void MeshPreprocessing::morphRBF2D()
{

}

void MeshPreprocessing::morphRBF3D(const std::string &name, const RBFTargetConfiguration &configuration)
{
	start("apply morphing '" + name + "'");

//	std::vector<double> pointsInMorphing;
//	std::vector<double> displacements;
//
//	std::vector<double> localPointsInMorphing;
//	std::vector<double> localDisplacements;
//
//	std::vector<double> w_matrix;
//	std::vector<double> q_matrix;
//
//	std::map<eslocal, std::pair<bool, std::vector<double>>> displacementForPoint;
//
//	// projde regiony a napocita body v moprhingu a posunuti
//	// localPointsInMorphing (x, y, z)/ (x, y)
//	// localDisplacements    (dx, dy, dy)/ (dx, dy)
//	for(auto region = configuration.targets.begin(); region!=configuration.targets.end();region++) {
//		ESINFO(OVERVIEW)<<"Processing region: "<<region->first;
//		Region* r = this->region(region->first);
//
//		Expression t_x(region->second.translation_x, {"x", "y", "z"});
//		Expression t_y(region->second.translation_y, {"x", "y", "z"});
//		Expression t_z(region->second.translation_z, {"x", "y", "z"});
//
//		for(auto element = r->elements().begin(); element!=r->elements().end(); element++){
//			//std::cout<<"Rank: "<<environment->MPIrank<<" - "<<**element<<" - nodes: "<<(*element)->nodes()<<std::endl;
//			for(size_t i = 0; i<(*element)->nodes();i++) {
//				eslocal node = (*element)->node(i);
//				if (this->nodes()[node]->clusters().front()==environment->MPIrank) {
//
//					//TODO: Improve this solution
//					auto elementInDisplacement = displacementForPoint.find(node);
//					if (elementInDisplacement==displacementForPoint.end()) {
//						std::pair<bool, std::vector<double>> empty;
//						empty.first = false;
//						empty.second.resize(3);
//
//						displacementForPoint[node] = empty;
//
//						elementInDisplacement = displacementForPoint.find(node);
//					}
//					if (elementInDisplacement->second.first && !region->second.overriding) {
//						//point overridden, but region is not overriding
//						continue;
//					}
//					elementInDisplacement->second.first=region->second.overriding;
//
//					const Point& p = this->coordinates()[node];
//
//					if (region->second.transformation == MORPHING_TRANSFORMATION::FIXED) {
//						for(int dim=0;dim<dimension;dim++) {
//							elementInDisplacement->second.second[dim]=0;
//						}
//					}else if (region->second.transformation == MORPHING_TRANSFORMATION::TRANSLATION) {
//						std::vector<double> point(3);
//						point[0] = p.x;
//						point[1] = p.y;
//						if (dimension==3){
//							point[2] = p.z;
//						}
//						elementInDisplacement->second.second[0]=(t_x.evaluate(point)-p.x);
//						elementInDisplacement->second.second[1]=(t_y.evaluate(point)-p.y);
//						if (dimension==3) {
//							elementInDisplacement->second.second[2]=(t_z.evaluate(point)-p.z);
//						}
//					}
//				}
//			}
//		}
//	}
//
//	for(auto it = displacementForPoint.begin(); it!=displacementForPoint.end(); it++) {
//		const Point& p = this->coordinates()[it->first];
//		localPointsInMorphing.push_back(p.x);
//		localPointsInMorphing.push_back(p.y);
//		if (dimension==3) {
//			localPointsInMorphing.push_back(p.z);
//		}
//		for(int dim=0;dim<dimension;dim++) {
//			localDisplacements.push_back( it->second.second[dim]);
//		}
//	}
//	// end
//	// localPointsInMorphing (x, y, z)/ (x, y)
//	// localDisplacements    (dx, dy, dy)/ (dx, dy)
//
//	displacementForPoint.clear();
//
//
//	// reseni se pocita na rank=0 ...
//	//std::cout<<"Process: "<<environment->MPIrank<<" sends "<<localPointsInMorphing.size()<<" coordinates\n";
//	Communication::gatherUnknownSize(localPointsInMorphing, pointsInMorphing);
//	//std::cout<<"Process: "<<environment->MPIrank<<" sends "<<localDisplacements.size()<<" displacements\n";
//	Communication::gatherUnknownSize(localDisplacements, displacements);
//
//	localPointsInMorphing.clear();
//	localDisplacements.clear();
//
//	// reseni se pocita na rank=0 -> tady se to provede
//	if (environment->MPIrank==0) {
//
//		int rowsFromCoordinates = displacements.size()/dimension;
//
//		DenseMatrix rhs(rowsFromCoordinates+dimension+1, dimension);
//		std::vector<double>::iterator it = displacements.begin();
//		for(int i=0;i<rowsFromCoordinates;i++) {
//			for(int k = 0;k<dimension;k++) {
//				rhs(i, k) = *it;
//				it++;
//			}
//		}
//
//		displacements.clear();
//
//		//std::cout<<rhs;
//
//		DenseMatrix M(rowsFromCoordinates + dimension+1, rowsFromCoordinates+dimension+1);
//
//		for (int i=0;i<rowsFromCoordinates;i++) {
//			for(int k = 0;k<dimension;k++) {
//				M(rowsFromCoordinates+k,i) = M(i,rowsFromCoordinates+k) = pointsInMorphing[i*dimension+k];
//			}
//			M(rowsFromCoordinates+dimension, i  ) =	M(i, rowsFromCoordinates+dimension) = 1;
//
//			for(int j=0;j<i;j++) {
//				double tmp = 0;
//				for(int k = 0;k<dimension;k++) {
//					double dist = pointsInMorphing[i*dimension+k] - pointsInMorphing[j*dimension+k];
//					tmp += dist*dist;
//				}
//				std::vector<double> function_data;
//				function_data.push_back(sqrt(tmp));
//				// TODO:
//				// M(i,j) = M(j,i) = configuration.function.evaluator->evaluate();
//			}
//		}
//
//		if (configuration.solver == MORPHING_RBF_SOLVER::DIRECT) {
//
//		} else if (configuration.solver == MORPHING_RBF_SOLVER::ITERATIVE) {
//
//			SparseCSRMatrix<MKL_INT> MSparse(M);
//			rhs.transpose();
//
//			DenseMatrix wq(rhs.rows(), rhs.columns());
//			for(int i=0;i<rhs.rows();i++) {
//				MSparse.gmresSolve(&rhs(i,0), &wq(i,0), configuration.solver_precision, 600);
//			}
//			//std::cout<<wq;
//
//			//TODO: optimize
//			for(int i=0;i<wq.columns();i++) {
//				for(int j=0;j<wq.rows(); j++) {
//					if (i<rowsFromCoordinates) {
//						w_matrix.push_back(wq(j,i));
//					}else {
//						q_matrix.push_back(wq(j,i));
//					}
//				}
//			}
//		}
//	}
//
//	// poslu vysledky na puvodni nody -> vsichni musi mit vsechno
//	Communication::broadcastUnknownSize(pointsInMorphing);
//	Communication::broadcastUnknownSize(w_matrix);
//	Communication::broadcastUnknownSize(q_matrix);
//
//	// zaloha coordinatu
//	_originalCoordinates = coordinates()._points;
//
//	int pointsCount = pointsInMorphing.size()/dimension;
//	for(int i=0;i<coordinates().clusterSize();i++) {
//		double displacement[dimension];
//
//		for(int k=0;k<dimension;k++) {
//			displacement[k]=0.0;
//		}
//		Point& p = (*_coordinates)[i];
//
//		std::vector<double> point;
//		point.push_back(p.x);
//		point.push_back(p.y);
//		if (dimension==3) point.push_back(p.z);
//		point.push_back(1);
//
//		for(int j=0;j<pointsCount;j++) {
//			double distance = 0.0;
//			for(int k=0;k<dimension;k++) {
//				double tmp = pointsInMorphing[j*dimension+k] - point[k];
//				distance+=tmp*tmp;
//			}
//			distance = sqrt(distance);
//
//			std::vector<double> function_data;
//			function_data.push_back(distance);
//			distance = function_r.evaluate(function_data);
//
//			for(int k=0;k<dimension;k++) {
//				displacement[k]+= distance * w_matrix[j*dimension + k];
//			}
//		}
//
//		for(int k=0;k<dimension;k++) {
//			for(int j=0;j<point.size();j++) {
//				displacement[k] += point[j]*q_matrix[j*dimension+k];
//			}
//		}
//
//		p.x += displacement[0];
//		p.y += displacement[1];
//		if (dimension==3) {
//			p.z += displacement[2];
//		}
//	}

	finish("apply morphing '" + name + "'");
}

