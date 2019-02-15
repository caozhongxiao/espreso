
#include "uniformnodescomposer.h"

#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/assembler/provider/provider.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/matrixtype.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "solver/generic/SparseMatrix.h"

#include <algorithm>
#include <numeric>


using namespace espreso;

void UniformNodesComposer::initDOFs()
{
	// ASSUME THAT SHARED NODES ARE SORTED IN THE SAME ORDER ON ALL PROCESSES

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<esint> doffset(threads);
	std::vector<std::vector<esint> > roffset(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = info::mesh->nodes->ranks->begin(t);
		esint dsize = 0;
		std::vector<esint> troffset(info::mesh->neighbours.size());

		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == info::mpi::rank) {
				dsize += _DOFs;
			} else {
				esint noffset = 0;
				while (info::mesh->neighbours[noffset] < ranks->front()) {
					++noffset;
				}
				++troffset[noffset];
			}
		}

		doffset[t] = dsize;
		roffset[t].swap(troffset);
	}

	esint goffset = Esutils::sizesToOffsets(doffset);
	Communication::exscan(goffset);
	Esutils::sizesToOffsets(roffset);

	std::vector<std::vector<std::vector<esint> > > sBuffer(threads, std::vector<std::vector<esint> >(info::mesh->neighbours.size()));
	std::vector<std::vector<esint> > rBuffer(info::mesh->neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = info::mesh->nodes->ranks->begin(t);
		esint toffset = doffset[t];
		std::vector<std::vector<esint> > tBuffer(info::mesh->neighbours.size());

		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == info::mpi::rank) {
				esint noffset = 0;
				for (auto r = ranks->begin() + 1; r != ranks->end(); ++r) {
					while (info::mesh->neighbours[noffset] < *r) {
						++noffset;
					}
					tBuffer[noffset].push_back(toffset + goffset);
				}
				toffset += _DOFs;
			}
		}
		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < sBuffer[t].size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	std::vector<std::vector<esint> > DOFDistribution(threads), DOFData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = info::mesh->nodes->ranks->begin(t);
		esint toffset = doffset[t];
		std::vector<esint> troffset = roffset[t];
		std::vector<esint> tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}

		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == info::mpi::rank) {
				for (int dof = 0; dof < _DOFs; ++dof) {
					tdata.push_back(goffset + toffset++);
				}
			} else {
				esint noffset = 0;
				while (info::mesh->neighbours[noffset] < ranks->front()) {
					++noffset;
				}
				for (int dof = 0; dof < _DOFs; ++dof) {
					tdata.push_back(rBuffer[noffset][troffset[noffset]] + dof);
				}
				++troffset[noffset];
			}
			tdist.push_back(tdata.size());
		}

		DOFDistribution[t].swap(tdist);
		DOFData[t].swap(tdata);
	}

	Esutils::threadDistributionToFullDistribution(DOFDistribution);

	_DOFMap = new serializededata<esint, esint>(DOFDistribution, DOFData);
}

void UniformNodesComposer::buildDirichlet()
{
	std::vector<std::vector<esint> > dIndices;
	_controler.dirichletIndices(dIndices);

	if (dIndices.size() != (size_t)_DOFs) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid number of DOFs per node in Dirichlet.";
	}

	if (_DOFs > 1) {
		#pragma omp parallel for
		for (int dof = 0; dof < _DOFs; ++dof) {
			for (size_t i = 0; i < dIndices[dof].size(); ++i) {
				dIndices[dof][i] = _DOFs * dIndices[dof][i] + dof;
			}
		}
	}
	for (int dof = 0; dof < _DOFs; ++dof) {
		_dirichletMap.insert(_dirichletMap.end(), dIndices[dof].begin(), dIndices[dof].end());
	}

	_dirichletPermutation.resize(_dirichletMap.size());
	std::iota(_dirichletPermutation.begin(), _dirichletPermutation.end(), 0);
	std::sort(_dirichletPermutation.begin(), _dirichletPermutation.end(), [&] (esint i, esint j) {
		return _dirichletMap[i] < _dirichletMap[j];
	});

	std::sort(_dirichletMap.begin(), _dirichletMap.end());
}

void UniformNodesComposer::buildPatterns()
{
	size_t threads = info::env::OMP_NUM_THREADS;
//	MatrixType mtype = _provider.getMatrixType();
	MatrixType mtype = MatrixType::REAL_UNSYMMETRIC; // TODO: set dirichlet in symmetric matrices

	_nDistribution = info::mesh->nodes->gatherUniqueNodeDistribution();
	for (size_t n = 0; n < _nDistribution.size(); ++n) {
		_nDistribution[n] *= _DOFs;
	}
	_foreignDOFs = std::lower_bound(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end(), _nDistribution[info::mpi::rank]) - _DOFMap->datatarray().begin();

	_tKOffsets.resize(threads);
	_tRHSOffsets.resize(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint tKsize = 0, tRHSsize = 0;
		for (auto e = info::mesh->elements->procNodes->begin(t); e != info::mesh->elements->procNodes->end(t); ++e) {
			tRHSsize += e->size() * _DOFs;
			tKsize += getMatrixSize(e->size() * _DOFs, mtype);
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				for (auto e = info::mesh->boundaryRegions[r]->procNodes->begin(t); e != info::mesh->boundaryRegions[r]->procNodes->end(t); ++e) {
					tRHSsize += e->size() * _DOFs;
					tKsize += getMatrixSize(e->size() * _DOFs, mtype);;
				}
			}

			_tKOffsets[t] = tKsize;
			_tRHSOffsets[t] = tRHSsize;
		}
	}

	_localKOffset = Esutils::sizesToOffsets(_tKOffsets);
	_localRHSOffset = Esutils::sizesToOffsets(_tRHSOffsets);

	std::vector<IJ> KPattern(_localKOffset);
	std::vector<esint> RHSPattern(_localRHSOffset);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		IJ *Koffset = KPattern.data() + _tKOffsets[t];
		esint *RHSoffset = RHSPattern.data() + _tRHSOffsets[t];

		auto insert = [&] (serializededata<esint, esint>::const_iterator &e) {
			esint *_RHS = RHSoffset;
			for (esint dof = 0; dof < _DOFs; ++dof) {
				for (auto n = e->begin(); n != e->end(); ++n, ++RHSoffset) {
					*RHSoffset = _DOFMap->datatarray()[*n * _DOFs + dof];
				}
			}
			insertKPattern(Koffset, _RHS, RHSoffset, mtype);
		};

		for (auto e = info::mesh->elements->procNodes->cbegin(t); e != info::mesh->elements->procNodes->cend(t); ++e) {
			insert(e);
			Koffset += getMatrixSize(e->size() * _DOFs, mtype);
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				for (auto e = info::mesh->boundaryRegions[r]->procNodes->cbegin(t); e != info::mesh->boundaryRegions[r]->procNodes->cend(t); ++e) {
					insert(e);
					Koffset += getMatrixSize(e->size() * _DOFs, mtype);
				}
			}
		}
	}

	std::vector<esint> pK(KPattern.size());
	std::iota(pK.begin(), pK.end(), 0);
	std::sort(pK.begin(), pK.end(), [&] (esint i, esint j) {
		return KPattern[i] < KPattern[j];
	});

	std::vector<esint> pRHS(RHSPattern.size());
	std::iota(pRHS.begin(), pRHS.end(), 0);
	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	std::vector<std::vector<IJ> > sKBuffer(info::mesh->neighbours.size()), rKBuffer(info::mesh->neighbours.size());
	std::vector<std::vector<esint> > sRHSBuffer(info::mesh->neighbours.size()), rRHSBuffer(info::mesh->neighbours.size());

	auto iK = pK.begin();
	auto iRHS = pRHS.begin();
	for (size_t n = 0; n < info::mesh->neighbours.size() && info::mesh->neighbours[n] < info::mpi::rank; ++n) {
		while (KPattern[*iK].row < _nDistribution[info::mesh->neighbours[n] + 1]) {
			if (iK == pK.begin() || KPattern[*iK] != KPattern[*(iK - 1)]) {
				sKBuffer[n].push_back(KPattern[*iK]);
			}
			++iK;
		}
		while (RHSPattern[*iRHS] < _nDistribution[info::mesh->neighbours[n] + 1]) {
			if (iRHS == pRHS.begin() || RHSPattern[*iRHS] != RHSPattern[*(iRHS - 1)]) {
				sRHSBuffer[n].push_back(RHSPattern[*iRHS]);
			}
			++iRHS;
		}
	}

	if (!Communication::receiveUpperUnknownSize(sKBuffer, rKBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange K pattern.";
	}
	if (!Communication::receiveUpperUnknownSize(sRHSBuffer, rRHSBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange RHS pattern.";
	}

	for (size_t i = 0; i < rKBuffer.size(); i++) {
		KPattern.insert(KPattern.end(), rKBuffer[i].begin(), rKBuffer[i].end());
	}
	for (size_t i = 0; i < rRHSBuffer.size(); i++) {
		RHSPattern.insert(RHSPattern.end(), rRHSBuffer[i].begin(), rRHSBuffer[i].end());
	}

	_nKSize.clear();
	_nRHSSize.clear();
	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		if (info::mesh->neighbours[n] < info::mpi::rank) {
			_nKSize.push_back(sKBuffer[n].size());
			_nRHSSize.push_back(sRHSBuffer[n].size());
		} else {
			_nKSize.push_back(rKBuffer[n].size());
			_nRHSSize.push_back(rRHSBuffer[n].size());
		}
	}

	size_t localK = pK.size(), localRHS = pRHS.size();
	pK.resize(KPattern.size());
	pRHS.resize(RHSPattern.size());
	std::iota(pK.begin() + localK, pK.end(), localK);
	std::iota(pRHS.begin() + localRHS, pRHS.end(), localRHS);
	std::sort(pK.begin(), pK.end(), [&] (esint i, esint j) {
		return KPattern[i] < KPattern[j];
	});
	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	_KPermutation.resize(KPattern.size());
	_RHSPermutation.resize(RHSPattern.size());
	data->K.resize(1);
	data->M.resize(1);
	data->f.resize(1);
	data->R.resize(1);
	data->primalSolution.resize(1);
	data->dualSolution.resize(1);
	data->K.front().haloRows = (info::mesh->nodes->size - info::mesh->nodes->uniqueSize) * _DOFs;
	data->K.front().rows = info::mesh->nodes->size * _DOFs;
	data->K.front().cols = info::mesh->nodes->uniqueTotalSize * _DOFs;
	data->f.front().resize(info::mesh->nodes->size * _DOFs);
	data->R.front().resize(info::mesh->nodes->size * _DOFs);
	data->primalSolution.front().resize(info::mesh->nodes->size * _DOFs);

	std::vector<esint> &ROW = data->K.front().CSR_I_row_indices;
	std::vector<esint> &COL = data->K.front().CSR_J_col_indices;
	std::vector<double> &VAL  = data->K.front().CSR_V_values;

	ROW.reserve(info::mesh->nodes->size * _DOFs + 1);
	ROW.push_back(1);
	COL.push_back(KPattern[pK.front()].column + 1);
	_KPermutation[pK.front()] = 0;
	_RHSPermutation[pRHS.front()] = 0;
	esint maxcol = 0, mincol = data->K.front().cols;
	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
			++nonzeros;
			COL.push_back(KPattern[pK[i]].column + 1);
			maxcol = std::max(maxcol, KPattern[pK[i]].column + 1);
			mincol = std::min(mincol, KPattern[pK[i]].column + 1);
			if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
				ROW.push_back(nonzeros + 1);
			}
		}
		_KPermutation[pK[i]] = nonzeros;
	}
	for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
			++nonzeros;
		}
		_RHSPermutation[pRHS[i]] = nonzeros;
	}
	ROW.push_back(COL.size() + 1);
	VAL.resize(COL.size());

	data->K.front().minCol = mincol;
	data->K.front().maxCol = maxcol;
	data->K.front().nnz = COL.size();
	data->K.front().mtype = _provider.getMatrixType();
	switch (data->K.front().mtype) {
	case MatrixType::REAL_UNSYMMETRIC: data->K.front().type = 'G'; break;
	default: data->K.front().type = 'S';
	}
	data->K.front().CSR_I_row_indices.swap(ROW);
	data->K.front().CSR_J_col_indices.swap(COL);
	data->M.front() = data->K.front();
	data->M.front().mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	data->M.front().type = 'S';
}

void UniformNodesComposer::synchronize(Matrices matrices)
{
	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	for (size_t n = 0; n < info::mesh->neighbours.size(); n++) {
		esint size = 0;
		size += (matrices & Matrices::K) ? _nKSize[n] : 0;
		size += (matrices & Matrices::M) ? _nKSize[n] : 0;
		size += (matrices & Matrices::R) ? _nRHSSize[n] : 0;
		size += (matrices & Matrices::f) ? _nRHSSize[n] : 0;
		if (info::mpi::rank < info::mesh->neighbours[n]) {
			rBuffer[n].resize(size);
		} else {
			sBuffer[n].reserve(size);
		}
	}

	auto nranks = info::mesh->nodes->ranks->begin();
	auto DOFs = _DOFMap->begin();
	if (matrices & Matrices::K) {
		for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
			esint r = 0;
			while (info::mesh->neighbours[r] < nranks->front()) {
				++r;
			}
			auto begin = data->K.front().CSR_I_row_indices[DOFs->begin() - _DOFMap->datatarray().begin()] - 1;
			auto end = data->K.front().CSR_I_row_indices[DOFs->end() - _DOFMap->datatarray().begin()] - 1;
			sBuffer[r].insert(sBuffer[r].end(), data->K.front().CSR_V_values.begin() + begin, data->K.front().CSR_V_values.begin() + end);
		}
		nranks = info::mesh->nodes->ranks->begin();
		DOFs = _DOFMap->begin();
	}
	if (matrices & Matrices::M) {
		for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
			esint r = 0;
			while (info::mesh->neighbours[r] < nranks->front()) {
				++r;
			}
			auto begin = data->M.front().CSR_I_row_indices[DOFs->begin() - _DOFMap->datatarray().begin()] - 1;
			auto end = data->M.front().CSR_I_row_indices[DOFs->end() - _DOFMap->datatarray().begin()] - 1;
			sBuffer[r].insert(sBuffer[r].end(), data->M.front().CSR_V_values.begin() + begin, data->M.front().CSR_V_values.begin() + end);
		}
		nranks = info::mesh->nodes->ranks->begin();
		DOFs = _DOFMap->begin();
	}
	if (matrices & Matrices::R) {
		for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
			esint r = 0;
			while (info::mesh->neighbours[r] < nranks->front()) {
				++r;
			}
			for (int dof = 0; dof < _DOFs; ++dof) {
				sBuffer[r].push_back(data->R.front()[n * _DOFs + dof]);
			}
		}
		nranks = info::mesh->nodes->ranks->begin();
		DOFs = _DOFMap->begin();
	}
	if (matrices & Matrices::f) {
		for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
			esint r = 0;
			while (info::mesh->neighbours[r] < nranks->front()) {
				++r;
			}
			for (int dof = 0; dof < _DOFs; ++dof) {
				sBuffer[r].push_back(data->f.front()[n * _DOFs + dof]);
			}
		}
	}


	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize assembled data.";
	}

	size_t KIndex = _localKOffset, RHSIndex = _localRHSOffset;
	for (size_t n = 0; n < info::mesh->neighbours.size(); n++) {
		if (info::mpi::rank < info::mesh->neighbours[n]) {
			size_t index = 0, k = KIndex, r = RHSIndex;
			if (matrices & Matrices::K) {
				for (esint i = 0; i < _nKSize[n]; ++k, ++index, ++i) {
					data->K.front().CSR_V_values[_KPermutation[k]] += rBuffer[n][index];
				}
			}
			k = KIndex;
			if (matrices & Matrices::M) {
				for (esint i = 0; i < _nKSize[n]; ++k, ++index, ++i) {
					data->M.front().CSR_V_values[_KPermutation[k]] += rBuffer[n][index];
				}
			}
			KIndex += _nKSize[n];
			if (matrices & Matrices::R) {
				for (esint i = 0; i < _nRHSSize[n]; ++r, ++index, ++i) {
					data->R.front()[_RHSPermutation[r]] += rBuffer[n][index];
				}
			}
			r = RHSIndex;
			if (matrices & Matrices::f) {
				for (esint i = 0; i < _nRHSSize[n]; ++r, ++index, ++i) {
					data->f.front()[_RHSPermutation[r]] += rBuffer[n][index];
				}
			}
			RHSIndex += _nRHSSize[n];
		}
	}
}

void UniformNodesComposer::gather(std::vector<double> &data)
{
	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	size_t RHSIndex = _localRHSOffset;
	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		if (info::mesh->neighbours[n] < info::mpi::rank) {
			rBuffer[n].resize(_nRHSSize[n]);
		} else {
			sBuffer[n].reserve(_nRHSSize[n]);
			for (esint i = 0; i < _nRHSSize[n]; ++i) {
				sBuffer[n].push_back(data[_RHSPermutation[RHSIndex++]]);
			}
		}
	}

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize assembled data.";
	}

	std::vector<esint> rIndices(info::mesh->neighbours.size());
	auto nranks = info::mesh->nodes->ranks->begin();
	auto DOFs = _DOFMap->begin();
	for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
		esint r = 0;
		while (info::mesh->neighbours[r] < nranks->front()) {
			++r;
		}
		for (int dof = 0; dof < _DOFs; ++dof) {
			data[n * _DOFs + dof] = rBuffer[r][rIndices[r]++];
		}
	}
}
