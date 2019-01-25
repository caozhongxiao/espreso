
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "physics/assembler/dataholder.h"
#include "physics/assembler/composer/composer.h"
#include "assembler.h"

#include "basis/logging/logging.h"
#include "basis/utilities/communication.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "solver/generic/SparseMatrix.h"
#include "linearsolver/linearsolver.h"

using namespace espreso;

void Assembler::callsolve(Matrices matrices)
{
	_solver->solve(matrices);
}

double Assembler::lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext)
{
//	double alpha = 1;
//	auto multiply = [] (const std::vector<std::vector<double> > &v1, const std::vector<std::vector<double> > &v2) {
//		double cmul = 0, gmul;
//
//		#pragma omp parallel for
//		for (size_t d = 0; d < v1.size(); d++) {
//			double dmul = 0;
//			for (size_t i = 0; i < v1[d].size(); i++) {
//				dmul += v1[d][i] * v2[d][i];
//			}
//			#pragma omp atomic
//			cmul += dmul;
//		}
//		MPI_Allreduce(&cmul, &gmul, 1, MPI_DOUBLE, MPI_SUM, info::mpi::MPICommunicator);
//		return gmul;
//	};
//
//	double a = 0, b = 1;
//	double fa = 0, fb = 0, fx = 0, faStart = 0;
//
//	std::vector<std::vector<double> > solution = deltaU;
//	std::vector<std::vector<double> > F_ext_r = F_ext;
//
//	for (size_t i = 0; i < 6; i++) {
//		sum(solution, 1, U, alpha, deltaU, "U = U + alpha * delta U (line search)");
//
//		solution.swap(info::data->primalSolution);
////			physics.updateMatrix(Matrices::R);
//		solution.swap(info::data->primalSolution);
//
//		if (i == 0) {
//			faStart = multiply(deltaU, info::data->f);
//			sum(F_ext_r, 1, F_ext, -1, info::data->R, "F_ext - R");
//
//			fb = multiply(deltaU, F_ext_r);
//			if ((faStart < 0 && fb < 0) || (faStart >= 0 && fb >= 0)) {
//				return alpha;
//			}
//			fa = faStart;
//		} else {
//			sum(F_ext_r, 1, F_ext, -1, info::data->R, "F_ext - R");
//			fx = multiply(deltaU, F_ext_r);
//			if (fa * fx < 0) {
//				b = alpha;
//				fb = fx;
//			} else if (fb * fx < 0) {
//				a = alpha;
//				fa = fx;
//			}
//
//			if (fabs(fx) <= 0.5 * faStart) {
//				alpha = a - fa * ((b - a ) / (fb - fa));
//				break;
//			}
//		}
//
//		alpha = a - fa * ((b - a ) / (fb - fa));
//	}
//
//	if (alpha < 0.1) {
//		alpha = 0.1;
//	}
//	if (alpha > .99) {
//		alpha = 1;
//	}
//
//	sum(solution, 0, U, alpha, deltaU, "delta U = alpha * delta U (line search)");
//	solution.swap(deltaU);
//	return alpha;
	return 0;
}





