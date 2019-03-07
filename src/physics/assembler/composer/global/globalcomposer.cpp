
#include "physics/assembler/dataholder.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/assembler/assembler.h"

#include "globalcomposer.h"

#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "solver/generic/SparseMatrix.h"

#include "wrappers/math/math.h"

#include <algorithm>

using namespace espreso;

void GlobalComposer::buildMVData()
{
	std::vector<esint> &ROWS = data->K.front().CSR_I_row_indices;
	std::vector<esint> &COLS = data->K.front().CSR_J_col_indices;

	std::vector<IJ> synchronization;

	_MVRows.push_back(1);
	for (esint r = data->K.front().haloRows; r < data->K.front().rows; ++r) {
		for (esint c = ROWS[r] - 1; c < ROWS[r + 1] - 1; ++c) {
			_MVCols.push_back(COLS[c] - data->K.front().minCol + 1);

			esint cindex = COLS[c] - 1;
			if (cindex < _nDistribution[info::mpi::rank] || _nDistribution[info::mpi::rank + 1] <= cindex) {
				synchronization.push_back({r, COLS[c] - 1});
			}
		}
		_MVRows.push_back(_MVCols.size() + 1);
	}
	_MVValuesOffset = ROWS[data->K.front().haloRows] - 1;
	_MVVec.resize(data->K.front().maxCol - data->K.front().minCol + 1);

	std::sort(synchronization.begin(), synchronization.end(), [] (const IJ &i, const IJ &j) {
		if (i.column == j.column) {
			return i.row < j.row;
		}
		return i.column < j.column;
	});

	_MVNeighbours = info::mesh->neighbours;
	auto sbegin = synchronization.begin();
	while (sbegin != synchronization.end()) {
		int neighbor = std::lower_bound(_nDistribution.begin(), _nDistribution.end(), sbegin->column + 1) - _nDistribution.begin() - 1;
		if (neighbor != info::mpi::rank) {
			_MVNeighbours.push_back(neighbor);
		}
		sbegin = std::lower_bound(sbegin, synchronization.end(), _nDistribution[neighbor + 1], [] (const IJ &index, esint value) {
			return index.column < value;
		});
	}
	utils::sortAndRemoveDuplicity(_MVNeighbours);

	_MVSend.resize(_MVNeighbours.size());
	_MVRecv.resize(_MVNeighbours.size());

	for (size_t n = 0, i = 0; n < _MVNeighbours.size(); n++) {
		while (i < synchronization.size() && synchronization[i].column < _nDistribution[_MVNeighbours[n] + 1]) {
			_MVSend[n].push_back(synchronization[i].row);
			_MVRecv[n].push_back(synchronization[i].column - data->K.front().minCol + 1);
			++i;
		}
	}

	for (size_t n = 0; n < _MVNeighbours.size(); n++) {
		utils::sortAndRemoveDuplicity(_MVSend[n]);
		utils::sortAndRemoveDuplicity(_MVRecv[n]);
	}
}

void GlobalComposer::assemble(Matrices matrices, const SolverParameters &parameters)
{
	if (!(matrices & (Matrices::K | Matrices::M | Matrices::R | Matrices::f))) {
		return;
	}

	size_t threads = info::env::OMP_NUM_THREADS;

//	MatrixType mtype = _provider.getMatrixType();
	MatrixType mtype = MatrixType::REAL_UNSYMMETRIC; // TODO: set dirichlet in symmetric matrix

	clearMatrices(matrices, 0);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t KIndex = _tKOffsets[t], RHSIndex = _tRHSOffsets[t];
		double KReduction = parameters.timeIntegrationConstantK, RHSReduction = parameters.internalForceReduction;
		Controller::InstanceFiller filler;

		switch (mtype) {
		case MatrixType::REAL_UNSYMMETRIC:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						#pragma omp atomic
						data->f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						#pragma omp atomic
						data->R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = 0; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							#pragma omp atomic
							data->K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows() && r / filler.Me.rows() == c / filler.Me.rows()) {
							#pragma omp atomic
							data->M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r % filler.Me.rows(), c % filler.Me.rows());
						}
						if ((matrices & Matrices::C) && filler.Ce.rows()) {
							#pragma omp atomic
							data->C[0].CSR_V_values[_KPermutation[KIndex]] += filler.Ce(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						#pragma omp atomic
						data->f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						#pragma omp atomic
						data->R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							#pragma omp atomic
							data->K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows() && r / filler.Me.rows() == c / filler.Me.rows()) {
							#pragma omp atomic
							data->M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r % filler.Me.rows(), c % filler.Me.rows());
						}
						if ((matrices & Matrices::C) && filler.Ce.rows()) {
							#pragma omp atomic
							data->C[0].CSR_V_values[_KPermutation[KIndex]] += filler.Ce(r, c);
						}
					}
				}
			}; break;
		}

		filler.begin = info::mesh->elements->distribution[t];
		filler.end = info::mesh->elements->distribution[t + 1];

		_controler.processElements(matrices, parameters, filler);

		KReduction = parameters.internalForceReduction;
		filler.Me.resize(0, 0);
		filler.Re.resize(0, 0);

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->distribution.size()) {
				filler.begin = info::mesh->boundaryRegions[r]->distribution[t];
				filler.end = info::mesh->boundaryRegions[r]->distribution[t + 1];
				_controler.processBoundary(matrices, parameters, r, filler);
			}
		}
	}
	synchronize(matrices);
}

struct __Dirichlet__ {
	double value;
	esint row, column;
};

void GlobalComposer::setDirichlet(Matrices matrices, double reduction, const std::vector<double> &subtraction)
{
	if (matrices & (Matrices::K | Matrices::M)) {
		_KDirichletValues.clear();
	}

	auto &ROW = data->K.front().CSR_I_row_indices;
	auto &COL = data->K.front().CSR_J_col_indices;
	auto &VAL = data->K.front().CSR_V_values;
	auto &RHS = data->f.front();

	std::vector<double> values(_dirichletMap.size());
	_controler.dirichletValues(values);

	for (size_t i = 0; i < values.size(); i++) {
		values[i] *= reduction;
	}

	if (subtraction.size()) {
		for (size_t i = 0; i < _dirichletMap.size(); i++) {
			values[_dirichletPermutation[i]] -= subtraction[_dirichletMap[i]];
		}
	}

	std::vector<__Dirichlet__> dirichlet;

	for (size_t i = 0, j = 0, v = 0; i < _dirichletMap.size(); i = j) {
		double value = 0;
		while (j < _dirichletMap.size() && _dirichletMap[j] == _dirichletMap[i]) {
			value += values[_dirichletPermutation[j++]];
		}
		value /= j - i;
		RHS[_dirichletMap[i]] = value;
		esint col = _DOFMap->datatarray()[_dirichletMap[i]] + 1;
		for (esint j = ROW[_dirichletMap[i]]; j < ROW[_dirichletMap[i] + 1]; j++) {
			if (COL[j - 1] == col) {
				VAL[j - 1] = 1;
			} else {
				VAL[j - 1] = 0;
				auto dit = std::lower_bound(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end(), COL[j - 1] - 1);
				if (
						(dit != _DOFMap->datatarray().end() && *dit == COL[j - 1] - 1) &&
						!std::binary_search(_dirichletMap.begin(), _dirichletMap.end(), dit - _DOFMap->datatarray().begin())
						) {
					esint r = dit - _DOFMap->datatarray().begin();
					for (esint c = ROW[r]; c < ROW[r + 1]; c++) {
						if (COL[c - 1] == col) {
							double val;
							if (v == _KDirichletValues.size()) {
								val = VAL[c - 1];
								_KDirichletValues.push_back(val);
							} else {
								val = _KDirichletValues[v];
							}
							++v;
							if (r < _foreignDOFs) {
								dirichlet.push_back({val * RHS[_dirichletMap[i]], *dit, COL[c - 1] - 1});
							}
							RHS[r] -= val * RHS[_dirichletMap[i]];
							VAL[c - 1] = 0;
						}
					}
				} else {
					// the column will be updated by a higher process
				}
			}
		}
	}

	std::sort(dirichlet.begin(), dirichlet.end(), [] (const __Dirichlet__ &i, const __Dirichlet__ &j) {
		if (i.row == j.row) {
			return i.column < j.column;
		}
		return i.row < j.row;
	});

	std::vector<std::vector<__Dirichlet__> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());
	for (size_t n = 0, i = 0; n < info::mesh->neighbours.size(); n++) {
		while (i < dirichlet.size() && dirichlet[i].row < _nDistribution[info::mesh->neighbours[n] + 1]) {
			sBuffer[n].push_back(dirichlet[i++]);
		}
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		eslog::error("ESPRESO internal error: synchronize dirichlet data.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbours.size(); n++) {
		for (size_t i = 0; i < rBuffer[n].size(); i++) {
			if (!std::binary_search(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end(), rBuffer[n][i].column)) {
				esint r = std::lower_bound(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end(), rBuffer[n][i].row) - _DOFMap->datatarray().begin();
				RHS[r] -= rBuffer[n][i].value;
				esint c = std::lower_bound(COL.begin() + ROW[r] - 1, COL.begin() + ROW[r + 1] - 1, rBuffer[n][i].column + 1) - COL.begin();
				VAL[c] = 0;
			}
		}
	}
}

void GlobalComposer::fillSolution()
{
	memcpy(
			_controler.solution()->data.data() + data->K.front().haloRows,
			data->primalSolution.front().data() + data->K.front().haloRows,
			sizeof(double) * (data->K.front().rows - data->K.front().haloRows));

	gather(_controler.solution()->data);
}

void GlobalComposer::KplusAlfaM(double alfa)
{
	data->K[0].MatAddInPlace(data->M[0], 'N', alfa);
}

void GlobalComposer::alfaKplusBetaM(double alfa, double beta)
{
	data->K[0].MatScale(alfa);
	KplusAlfaM(beta);
}

void GlobalComposer::apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x)
{
	esint rows = data->K.front().rows - data->K.front().haloRows;
	esint cols = data->K.front().maxCol - data->K.front().minCol + 1;
	std::vector<double> &VALS = matrices.front().CSR_V_values;

	std::vector<std::vector<double> > sBuffer(_MVNeighbours.size()), rBuffer(_MVNeighbours.size());
	for (size_t n = 0; n < _MVNeighbours.size(); ++n) {
		rBuffer[n].resize(_MVRecv[n].size());
		sBuffer[n].reserve(_MVSend[n].size());
		for (size_t i = 0; i < _MVSend[n].size(); ++i) {
			sBuffer[n].push_back(x[_MVSend[n][i]]);
		}
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, _MVNeighbours)) {
		eslog::error("ESPRESO internal error: exchange MV data.\n");
	}

	for (size_t n = 0; n < _MVNeighbours.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			_MVVec[_MVRecv[n][i]] = rBuffer[n][i];
		}
	}
	std::copy(x.begin() + _foreignDOFs, x.end(), _MVVec.begin() + _nDistribution[info::mpi::rank] - data->K.front().minCol + 1);
	for (size_t n = 0; n < _MVNeighbours.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			_MVVec[_MVRecv[n][i]] = rBuffer[n][i];
		}
	}

	result.resize(data->K.front().rows);
	MATH::CSRMatVecProduct(rows, cols, _MVRows.data(), _MVCols.data(), VALS.data() + _MVValuesOffset, _MVVec.data(), result.data() + data->K.front().haloRows);

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

