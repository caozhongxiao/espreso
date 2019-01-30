
#include "physics/assembler/dataholder.h"
#include "physics/assembler/controllers/controller.h"
#include "globalcomposer.h"

#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "basis/containers/serializededata.h"
#include "basis/logging/logging.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "mesh/store/nodestore.h"
#include "solver/generic/SparseMatrix.h"

#include "wrappers/math/math.h"

#include <algorithm>

using namespace espreso;

void GlobalComposer::KplusAlfaM(double alfa)
{
	data->K[0].MatAddInPlace(data->M[0], 'N', alfa);
}

void GlobalComposer::apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x)
{
	// halo rows cannot be used since MKL does not accept it
	esint rows = data->K.front().rows;// - data->K.front().haloRows;
	esint cols = data->K.front().maxCol - data->K.front().minCol + 1;
	std::vector<esint> &ROWS = matrices.front().CSR_I_row_indices;
	std::vector<double> &VALS = matrices.front().CSR_V_values;

	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());
	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		rBuffer[n].resize(_MVRecv[n].size());
		sBuffer[n].reserve(_MVSend[n].size());
		for (size_t i = 0; i < _MVSend[n].size(); ++i) {
			sBuffer[n].push_back(x[_MVSend[n][i]]);
		}
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange MV data.";
	}

	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			_MVVec[_MVRecv[n][i]] = rBuffer[n][i];
		}
	}
	std::copy(x.begin() + _foreignDOFs, x.end(), _MVVec.begin() + _nDistribution[info::mpi::rank] - data->K.front().minCol + 1);
	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			_MVVec[_MVRecv[n][i]] = rBuffer[n][i];
		}
	}

	result.resize(_MVVec.size());
	MATH::CSRMatVecProduct(rows, cols, ROWS.data(), _MVCols.data(), VALS.data(), _MVVec.data(), result.data());

	gather(result);
}

void GlobalComposer::enrichRHS(double alfa, NodeData* x)
{
	for (size_t i = 0; i < data->f[0].size(); ++i) {
		data->f[0][i] += alfa * x->data[i];
	}
}

void GlobalComposer::computeReactionForces()
{
	apply(data->solverK, data->dualSolution.front(), _controler.solution()->data);
	for (size_t i = 0; i < data->dualSolution.front().size(); i++) {
		data->dualSolution.front()[i] = data->solverF.front()[i] - data->dualSolution.front()[i];
	}
}

double GlobalComposer::residualNormNumerator()
{
	double square = 0;
	for (size_t i = _foreignDOFs; i < data->f.front().size(); i++) {
		square += (data->f.front()[i] - data->dualSolution.front()[i]) * (data->f.front()[i] - data->dualSolution.front()[i]);
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}

double GlobalComposer::residualNormDenominator()
{
	double square = 0;
	for (size_t i = _foreignDOFs, d = 0; i < data->origF.front().size(); i++) {
		while (d < _dirichletMap.size() && (size_t)_dirichletMap[d] < i) { d++; }
		if (d == _dirichletMap.size() || (size_t)_dirichletMap[d] != i) {
			square += data->origF.front()[i] * data->origF.front()[i];
		} else {
			square += data->R.front()[i] * data->R.front()[i];
		}
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}

