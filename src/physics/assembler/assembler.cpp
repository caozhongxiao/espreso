
#include "assembler.h"

#include "physics/dataholder.h"

#include "globals/run.h"
#include "basis/logging/logging.h"
#include "basis/utilities/communication.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "solver/generic/SparseMatrix.h"

using namespace espreso;

void Assembler::keepK()
{
	#pragma omp parallel for
	for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
		run::data->origK[d] = run::data->K[d];
	}
}

void Assembler::sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description)
{
	std::vector<size_t> prefix(x.size());
	for (size_t i = 0; i < x.size(); i++) {
		prefix[i] = x[i].size();
	}
	sum(z, a, x, b, y, prefix, description);
}

void Assembler::sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description)
{
	if (z.size() == 0) {
		z.resize(x.size());
	}
	#pragma omp parallel for
	for (size_t d = 0; d < x.size(); d++) {
		if (z[d].size() == 0) {
			z[d].resize(x[d].size());
		}
		if (x[d].size() != y[d].size() || z[d].size() != x[d].size()) {
			ESINFO(ERROR) << "ESPRESO internal error while " << description << ". Vectors have different dimension.";
		}
		for (size_t i = 0; i < x[d].size() && i < prefix[d]; i++) {
			z[d][i] = a * x[d][i] + b * y[d][i];
		}
	}
}


/// y = A * x
void Assembler::multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description)
{
	#pragma omp parallel for
	for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
		A[d].MatVec(x[d], y[d], 'N', 0, 0, 0);
	}
}

// v = x * y
double Assembler::multiply(std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, const std::string &description)
{
	double sum = 0;
	double psum = 0;
	#pragma omp parallel for reduction(+:psum)
	for (size_t d = 0; d < x.size(); d++) {
		if (x[d].size() != y[d].size()) {
			ESINFO(ERROR) << "ESPRESO internal error while " << description << ". Vectors have different dimension.";
		}
		for (size_t i = 0; i < x[d].size(); i++) {
			psum += x[d][i] * y[d][i];
		}
	}
	sum = psum;
	return sum;
}

double Assembler::sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction, const std::string &description)
{
	double result;
//	result = physics.sumSquares(data, restriction);
	return result;
}

void Assembler::addToDirichletInB1(double a, const std::vector<std::vector<double> > &x)
{
	#pragma omp parallel for
	for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
		for (size_t j = 0; j < run::data->B1[d].J_col_indices.size(); j++) {
			if (run::data->B1[d].I_row_indices[j] > (esint)run::data->block[DataHolder::CONSTRAINT::DIRICHLET]) {
				break;
			}
			run::data->B1c[d][j] += a * x[d][run::data->B1[d].J_col_indices[j] - 1];
		}
	}
}

double Assembler::maxAbsValue(const std::vector<std::vector<double> > &v, const std::string &description)
{
	double gmax;
	double max = 0;
	for (size_t p = 0; p < v.size(); p++) {
		max = std::max(max, std::fabs(*std::max_element(v[p].begin(), v[p].end(), [] (const double v1, const double v2) { return std::fabs(v1) < std::fabs(v2); })));
	}

	MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);
	return gmax;
}

double Assembler::lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext)
{
	double alpha = 1;
	auto multiply = [] (const std::vector<std::vector<double> > &v1, const std::vector<std::vector<double> > &v2) {
		double cmul = 0, gmul;

		#pragma omp parallel for
		for (size_t d = 0; d < v1.size(); d++) {
			double dmul = 0;
			for (size_t i = 0; i < v1[d].size(); i++) {
				dmul += v1[d][i] * v2[d][i];
			}
			#pragma omp atomic
			cmul += dmul;
		}
		MPI_Allreduce(&cmul, &gmul, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
		return gmul;
	};

	double a = 0, b = 1;
	double fa = 0, fb = 0, fx = 0, faStart = 0;

	std::vector<std::vector<double> > solution = deltaU;
	std::vector<std::vector<double> > F_ext_r = F_ext;

	for (size_t i = 0; i < 6; i++) {
		sum(solution, 1, U, alpha, deltaU, "U = U + alpha * delta U (line search)");

		solution.swap(run::data->primalSolution);
//			physics.updateMatrix(Matrices::R);
		solution.swap(run::data->primalSolution);

		if (i == 0) {
			faStart = multiply(deltaU, run::data->f);
			sum(F_ext_r, 1, F_ext, -1, run::data->R, "F_ext - R");

			fb = multiply(deltaU, F_ext_r);
			if ((faStart < 0 && fb < 0) || (faStart >= 0 && fb >= 0)) {
				return alpha;
			}
			fa = faStart;
		} else {
			sum(F_ext_r, 1, F_ext, -1, run::data->R, "F_ext - R");
			fx = multiply(deltaU, F_ext_r);
			if (fa * fx < 0) {
				b = alpha;
				fb = fx;
			} else if (fb * fx < 0) {
				a = alpha;
				fa = fx;
			}

			if (fabs(fx) <= 0.5 * faStart) {
				alpha = a - fa * ((b - a ) / (fb - fa));
				break;
			}
		}

		alpha = a - fa * ((b - a ) / (fb - fa));
	}

	if (alpha < 0.1) {
		alpha = 0.1;
	}
	if (alpha > .99) {
		alpha = 1;
	}

	sum(solution, 0, U, alpha, deltaU, "delta U = alpha * delta U (line search)");
	solution.swap(deltaU);
	return alpha;
}





