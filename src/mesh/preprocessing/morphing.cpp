
#include "meshpreprocessing.h"

using namespace espreso;

void MeshPreprocessing::morph2D()
{

}

void MeshPreprocessing::morph3D()
{
//	void Mesh::morph(const MeshMorphing &morphing)
//	{
//		if (morphing.type == MORPHING_TYPE::RBF) {
//			//TODO: Repair dimension checking
//			int dimension = Point::dimension();
//			dimension=2;
//
//
//			for(auto morphingElement = morphing.rbf.begin(); morphingElement!=morphing.rbf.end();morphingElement++) {
//				ESINFO(OVERVIEW)<<"Processing morphing: "<<morphingElement->first;
//
//				Expression function_r(morphingElement->second.function, {"R"});
//
//				std::vector<double> pointsInMorphing;
//				std::vector<double> displacements;
//
//				std::vector<double> localPointsInMorphing;
//				std::vector<double> localDisplacements;
//
//				std::vector<double> w_matrix;
//				std::vector<double> q_matrix;
//
//				std::map<eslocal, std::vector<double>> displacementForPoint;
//				std::map<eslocal, int> displacementSetTimes;
//
//				for(auto region = morphingElement->second.targets.begin(); region!=morphingElement->second.targets.end();region++) {
//					ESINFO(OVERVIEW)<<"Processing region: "<<region->first;
//					Region* r = this->region(region->first);
//
//					Expression t_x(region->second.translation_x, {"x", "y", "z"});
//					Expression t_y(region->second.translation_y, {"x", "y", "z"});
//					Expression t_z(region->second.translation_z, {"x", "y", "z"});
//
//					for(auto element = r->elements().begin(); element!=r->elements().end(); element++){
//						//std::cout<<"Rank: "<<environment->MPIrank<<" - "<<**element<<" - nodes: "<<(*element)->nodes()<<std::endl;
//						for(size_t i = 0; i<(*element)->nodes();i++) {
//							eslocal node = (*element)->node(i);
//							if (this->nodes()[node]->clusters().front()==environment->MPIrank) {
//								const Point& p = this->coordinates()[node];
//
//								/*
//								Element *nodeElement = this->nodes()[node];
//								//3d
//								nodeElement->parentFaces();
//								//2d
//								nodeElement->parentFaces();
//								 */
//								//TODO: Improve this solution
//								auto elementInDisplacement = displacementForPoint.find(node);
//								if (elementInDisplacement==displacementForPoint.end() || region->second.overriding) {
//									std::vector<double> values(3);
//									for(int dim=0;dim<dimension;dim++) {
//										values.push_back(0);
//									}
//									displacementForPoint[node] = values;
//									displacementSetTimes[node] = 0;
//
//									elementInDisplacement = displacementForPoint.find(node);
//								}
//
//								if (region->second.transformation == MORPHING_TRANSFORMATION::FIXED) {
//									for(int dim=0;dim<dimension;dim++) {
//										displacementForPoint[node][dim]+=0;
//										displacementSetTimes[node]+= 1;
//									}
//								}else if (region->second.transformation == MORPHING_TRANSFORMATION::TRANSLATION) {
//									std::vector<double> point(3);
//									point[0] = p.x;
//									point[1] = p.y;
//									if (dimension==3){
//										point[2] = p.z;
//									}
//									displacementForPoint[node][0]+=(t_x.evaluate(point)-p.x);
//									displacementForPoint[node][1]+=(t_y.evaluate(point)-p.y);
//									if (dimension==3) {
//										displacementForPoint[node][1]+=(t_z.evaluate(point)-p.z);
//									}
//									displacementSetTimes[node] += 1;
//								}
//							}
//						}
//					}
//				}
//
//				/*printf("All displacements for points\n");
//				for(auto it = displacementForPoint.begin(); it!=displacementForPoint.end(); it++) {
//					printf("%d -> %f %f -> %d\n", it->first, it->second[0], it->second[1], displacementSetTimes[it->first] );
//				}*/
//
//				for(auto it = displacementForPoint.begin(); it!=displacementForPoint.end(); it++) {
//					const Point& p = this->coordinates()[it->first];
//					localPointsInMorphing.push_back(p.x);
//					localPointsInMorphing.push_back(p.y);
//					if (dimension==3) {
//						localPointsInMorphing.push_back(p.z);
//					}
//					int divideBy = displacementSetTimes[it->first];
//					for(int dim=0;dim<dimension;dim++) {
//						localDisplacements.push_back( it->second[dim]);
//					}
//				}
//
//				//std::cout<<"Process: "<<environment->MPIrank<<" sends "<<localPointsInMorphing.size()<<" coordinates\n";
//				Communication::gatherUnknownSize(localPointsInMorphing, pointsInMorphing);
//				//std::cout<<"Process: "<<environment->MPIrank<<" sends "<<localDisplacements.size()<<" displacements\n";
//				Communication::gatherUnknownSize(localDisplacements, displacements);
//
//				if (environment->MPIrank==0) {
//					/*std::cout<<"Global Indexes: "<<globalIndexes.size()<<std::endl;
//					for(auto x = globalIndexes.begin(); x!=globalIndexes.end(); x++) {
//						std::cout<<" "<<(*x);
//					}
//					std::cout<<"\n";
//					 */
//					std::cout<<"Coordinates: "<<pointsInMorphing.size()<<std::endl;
//					for(auto x = pointsInMorphing.begin(); x!=pointsInMorphing.end(); x++) {
//						std::cout<<" "<<(*x);
//					}
//					std::cout<<"\n";
//
//					std::cout<<"Displacement: "<<displacements.size()<<std::endl;
//					for(auto x = displacements.begin(); x!=displacements.end(); x++) {
//						std::cout<<" "<<(*x);
//					}
//					std::cout<<"\n";
//
//					int rowsFromCoordinates = displacements.size()/dimension;
//
//					std::cout<<"Number rows: "<<rowsFromCoordinates<<"\n";
//
//					DenseMatrix rhs(rowsFromCoordinates+dimension+1, dimension);
//					std::vector<double>::iterator it = displacements.begin();
//					for(int i=0;i<rowsFromCoordinates;i++) {
//						for(int k = 0;k<dimension;k++) {
//							rhs(i, k) = *it;
//							it++;
//						}
//					}
//
//					displacements.clear();
//
//					std::cout<<rhs;
//
//					DenseMatrix M(rowsFromCoordinates + dimension+1, rowsFromCoordinates+dimension+1);
//
//					for (int i=0;i<rowsFromCoordinates;i++) {
//						for(int k = 0;k<dimension;k++) {
//							M(rowsFromCoordinates+k,i) = M(i,rowsFromCoordinates+k) = pointsInMorphing[i*dimension+k];
//						}
//						M(rowsFromCoordinates+dimension, i  ) =	M(i, rowsFromCoordinates+dimension) = 1;
//
//						for(int j=0;j<i;j++) {
//							//std::cout<<i<<","<<j<<" ";
//							double tmp = 0;
//							for(int k = 0;k<dimension;k++) {
//								//std::cout<<pointsInMorphing[i*dimension+k]<<"-"<< pointsInMorphing[j*dimension+k]<<" ";
//								double dist = pointsInMorphing[i*dimension+k] - pointsInMorphing[j*dimension+k];
//								tmp += dist*dist;
//							}
//							//std::cout<<" = "<<sqrt(tmp)<<"\n";
//							std::vector<double> function_data;
//							function_data.push_back(sqrt(tmp));
//							M(i,j) = M(j,i) = function_r.evaluate(function_data);
//						}
//					}
//
//					std::cout<<M;
//
//					if (morphingElement->second.solver == MORPHING_RBF_SOLVER::DIRECT) {
//
//	/*			DenseMatrix MT;
//				MT = M;
//				MT.transpose();
//				rhs = MT*rhs;
//
//				M=MT*M;
//
//				std::cout<<M;
//				std::cout<<rhs;
//
//				std::vector<lapack_int> ipiv(M.rows());
//				lapack_int info;
//
//				info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', M.rows(), M.values(), M.rows());
//				std::cout<<"Result  "<<info<<"\n";
//
//				rhs.transpose();
//
//				info = LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', M.rows(),1, M.values(), M.rows(), rhs.values(), 1);
//				std::cout<<"Result  "<<info<<"\n";
//				std::cout<<rhs;
//				for(int i=0;i<rhs.columns();i++) {
//					for(int j=0;j<rhs.rows(); j++) {
//						if (i<rowsFromCoordinates) {
//							w_matrix.push_back(rhs(j,i));
//						}else {
//							q_matrix.push_back(rhs(j,i));
//						}
//					}
//				}
//
//	*/
//
//	/*			std::vector<MKL_INT> ipiv(M.rows());
//				MKL_INT info;
//				//char U = 'U';
//				//std::vector<double> work(M.rows());
//				//int rows = M.rows();
//				//int workl = M.rows(); //-1 to compute
//				info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', M.rows(), M.values(), M.rows(), &ipiv[0]);
//				//dsytrf( &U, &rows , M.values() , &rows, &ipiv[0], &work[0], &workl, &info );
//				std::cout<<"Result  "<<info<<"\n";
//				//printf("work  %f %f\n", work[0], work[1]);
//				rhs.transpose();
//
//				//dsytrs( &U, &rows, &dimension, M.values(), &rows, &ipiv[0], rhs.values(), &rows, &info );
//				info = LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'U', M.rows(),1, M.values(), M.rows(), &ipiv[0], rhs.values(), 1);
//				std::cout<<"Result  "<<info<<"\n";
//				std::cout<<rhs;
//				for(int i=0;i<rhs.columns();i++) {
//					for(int j=0;j<rhs.rows(); j++) {
//						if (i<rowsFromCoordinates) {
//							w_matrix.push_back(rhs(j,i));
//						}else {
//							q_matrix.push_back(rhs(j,i));
//						}
//					}
//				}
//				*/
//					/*char U = 'U';
//					eslocal info;
//					std::vector<eslocal>  m_ipiv;
//					dsptrf( &U, &m_cols, &m_dense_values[0], &m_ipiv[0] , &info );*/
//
//				/*espreso::SparseMatrix A;
//				A.type = 'S';
//				A.mtype = espreso::MatrixType::REAL_SYMMETRIC_INDEFINITE;
//				A.rows = numVertsCage1 + 4;
//				// number of rows
//				A.cols = numVertsCage1 + 4;
//				// number of columns
//				// fill A, n_rhs, rhs
//
//				std::vector<double> sol;
//
//				sol.resize(rhs.size());
//				A.ConvertCSRToDense(1);
//
//				espreso::DenseSolverMKL denseSolverMKL;
//				denseSolverMKL.ImportMatrix(A);
//
//				std::stringstream ss;
//				ss << "Solve";
//				denseSolverMKL.Factorization(ss.str());
//
//				denseSolverMKL.Solve(rhs, sol, n_rhs);*/
//
//					} else if (morphingElement->second.solver == MORPHING_RBF_SOLVER::ITERATIVE) {
//
//						SparseCSRMatrix<MKL_INT> MSparse(M);
//						rhs.transpose();
//
//						DenseMatrix wq(rhs.rows(), rhs.columns());
//						for(int i=0;i<rhs.rows();i++) {
//							MSparse.gmresSolve(&rhs(i,0), &wq(i,0), morphingElement->second.solver_precision, 600);
//						}
//						std::cout<<wq;
//
//						//TODO: optimize
//						for(int i=0;i<wq.columns();i++) {
//							for(int j=0;j<wq.rows(); j++) {
//								if (i<rowsFromCoordinates) {
//									w_matrix.push_back(wq(j,i));
//								}else {
//									q_matrix.push_back(wq(j,i));
//								}
//							}
//						}
//					}
//				}
//
//				Communication::broadcastUnknownSize(pointsInMorphing);
//				Communication::broadcastUnknownSize(w_matrix);
//				Communication::broadcastUnknownSize(q_matrix);
//
//				std::cout<<"Rank: "<<environment->MPIrank<<" points: "<<pointsInMorphing.size()<<std::endl;
//				for(auto x = pointsInMorphing.begin(); x!=pointsInMorphing.end(); x++) {
//					std::cout<<" "<<(*x);
//				}
//				std::cout<<"\n";
//
//				std::cout<<"Rank: "<<environment->MPIrank<<" w: "<<w_matrix.size()<<std::endl;
//				for(auto x = w_matrix.begin(); x!=w_matrix.end(); x++) {
//					std::cout<<" "<<(*x);
//				}
//				std::cout<<"\n";
//
//				std::cout<<"Rank: "<<environment->MPIrank<<" q: "<<q_matrix.size()<<std::endl;
//				for(auto x = q_matrix.begin(); x!=q_matrix.end(); x++) {
//					std::cout<<" "<<(*x);
//				}
//				std::cout<<"\n";
//
//				_originalCoordinates = coordinates()._points;
//
//				int pointsCount = pointsInMorphing.size()/dimension;
//				for(int i=0;i<coordinates().clusterSize();i++) {
//					double displacement[dimension];
//
//					for(int k=0;k<dimension;k++) {
//						displacement[k]=0.0;
//					}
//					Point& p = (*_coordinates)[i];
//
//					std::vector<double> point;
//					point.push_back(p.x);
//					point.push_back(p.y);
//					if (dimension==3) point.push_back(p.z);
//					point.push_back(1);
//
//					for(int j=0;j<pointsCount;j++) {
//						double distance = 0.0;
//						for(int k=0;k<dimension;k++) {
//							double tmp = pointsInMorphing[j*dimension+k] - point[k];
//							distance+=tmp*tmp;
//						}
//						distance = sqrt(distance);
//
//						std::vector<double> function_data;
//						function_data.push_back(distance);
//						distance = function_r.evaluate(function_data);
//
//						for(int k=0;k<dimension;k++) {
//							displacement[k]+= distance * w_matrix[j*dimension + k];
//						}
//					}
//
//					for(int k=0;k<dimension;k++) {
//						for(int j=0;j<point.size();j++) {
//							//printf(" %d %d - %f * %f \n",k,j,point[j],q_matrix[j*dimension+k]);
//							displacement[k] += point[j]*q_matrix[j*dimension+k];
//						}
//					}
//
//					p.x += displacement[0];
//					p.y += displacement[1];
//					if (dimension==3) {
//						p.z += displacement[2];
//						//printf("%f %f %f\n", displacement[0], displacement[1], displacement[2]);
//					}else {
//						//printf("%f %f\n", displacement[0], displacement[1]);
//					}
//				}
//
//				for(int i=0;i<coordinates().clusterSize();i++) {
//					Point& p1 = (*_coordinates)[i];
//					Point& p2 = (_originalCoordinates)[i];
//					printf("%f %f %f - %f %f %f\n", p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
//				}
//			}
//		}
//	}
}

