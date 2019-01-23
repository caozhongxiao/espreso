
#include "physics/assembler/dataholder.h"
#include "globalcomposer.h"

#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "basis/logging/logging.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "mesh/store/nodestore.h"
#include "solver/generic/SparseMatrix.h"

#include "wrappers/math/math.h"

#include <algorithm>

using namespace espreso;

NodeData* GlobalComposer::RHS()
{
	return NULL;
}

void GlobalComposer::KplusAlfaM(double alfa)
{
	data->K[0].MatAddInPlace(data->M[0], 'N', alfa);
}

void GlobalComposer::apply(std::vector<SparseMatrix> &matrices, NodeData *result, NodeData *x)
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
			sBuffer[n].push_back(x->data[_MVSend[n][i]]);
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
	std::copy(x->data.begin() + _foreignDOFs, x->data.end(), _MVVec.begin() + _nDistribution[info::mpi::rank] - data->K.front().minCol + 1);
	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			_MVVec[_MVRecv[n][i]] = rBuffer[n][i];
		}
	}

	MATH::CSRMatVecProduct(rows, cols, ROWS.data(), _MVCols.data(), VALS.data(), _MVVec.data(), result->data.data());

	gather(result);
}

void GlobalComposer::enrichRHS(double alfa, NodeData* a)
{
	for (size_t i = 0; i < data->f[0].size(); ++i) {
		data->f[0][i] += alfa * a->data[i];
	}
}

void GlobalComposer::RHSMinusR()
{

}

void GlobalComposer::DirichletMinusRHS()
{

}

double GlobalComposer::residualNorm()
{
	return std::sqrt(0);
}

