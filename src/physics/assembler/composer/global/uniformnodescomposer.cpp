
#include "uniformnodescomposer.h"

#include "../../controlers/controler.h"

#include "../../../instance.h"
#include "../../../step.h"
#include "../../../../basis/containers/serializededata.h"
#include "../../../../basis/matrices/matrixtype.h"
#include "../../../../basis/utilities/communication.h"
#include "../../../../basis/utilities/utils.h"
#include "../../../../config/ecf/environment.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"
#include "../../../../mesh/store/nodestore.h"
#include "../../../../mesh/store/boundaryregionstore.h"

#include "../../../../solver/generic/SparseMatrix.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void UniformNodesComposer::initDOFs()
{
	// ASSUME THAT SHARED NODES ARE SORTED IN THE SAME ORDER ON ALL PROCESSES

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> doffset(threads);
	std::vector<std::vector<eslocal> > roffset(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = _mesh.nodes->ranks->begin(t);
		eslocal dsize = 0;
		std::vector<eslocal> troffset(_mesh.neighbours.size());

		for (eslocal n = _mesh.nodes->distribution[t]; n < _mesh.nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				dsize += _DOFs;
			} else {
				eslocal noffset = 0;
				while (_mesh.neighbours[noffset] < ranks->front()) {
					++noffset;
				}
				++troffset[noffset];
			}
		}

		doffset[t] = dsize;
		roffset[t].swap(troffset);
	}

	eslocal goffset = Esutils::sizesToOffsets(doffset);
	Communication::exscan(goffset);
	Esutils::sizesToOffsets(roffset);

	std::vector<std::vector<std::vector<eslocal> > > sBuffer(threads, std::vector<std::vector<eslocal> >(_mesh.neighbours.size()));
	std::vector<std::vector<eslocal> > rBuffer(_mesh.neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = _mesh.nodes->ranks->begin(t);
		eslocal toffset = doffset[t];
		std::vector<std::vector<eslocal> > tBuffer(_mesh.neighbours.size());

		for (eslocal n = _mesh.nodes->distribution[t]; n < _mesh.nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				eslocal noffset = 0;
				for (auto r = ranks->begin() + 1; r != ranks->end(); ++r) {
					while (_mesh.neighbours[noffset] < *r) {
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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh.neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	std::vector<std::vector<eslocal> > DOFDistribution(threads), DOFData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = _mesh.nodes->ranks->begin(t);
		eslocal toffset = doffset[t];
		std::vector<eslocal> troffset = roffset[t];
		std::vector<eslocal> tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}

		for (eslocal n = _mesh.nodes->distribution[t]; n < _mesh.nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				for (int dof = 0; dof < _DOFs; ++dof) {
					tdata.push_back(goffset + toffset++);
				}
			} else {
				eslocal noffset = 0;
				while (_mesh.neighbours[noffset] < ranks->front()) {
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

	_DOFMap = new serializededata<eslocal, eslocal>(DOFDistribution, DOFData);

//	std::cout << environment->MPIrank << ": " << *_DOFMap << "\n";
}

void UniformNodesComposer::initDirichlet()
{
	std::vector<std::vector<eslocal> > dIndices;
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
	std::sort(_dirichletPermutation.begin(), _dirichletPermutation.end(), [&] (eslocal i, eslocal j) {
		return _dirichletMap[i] < _dirichletMap[j];
	});

	std::sort(_dirichletMap.begin(), _dirichletMap.end());
}

void UniformNodesComposer::buildPatterns()
{
	size_t threads = environment->OMP_NUM_THREADS;
	// MatrixType mtype = _controler.getMatrixType();
	MatrixType mtype = MatrixType::REAL_UNSYMMETRIC; // HYPRE not support symmetric systems

	std::vector<std::vector<eslocal> > RHSsize(threads), Ksize(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		eslocal tKsize = 0, tRHSsize = 0;
		for (auto e = _mesh.elements->procNodes->begin(t); e != _mesh.elements->procNodes->end(t); ++e) {
			tRHSsize += e->size() * _DOFs;
			tKsize += getMatrixSize(e->size() * _DOFs, mtype);
		}

		Ksize[t].push_back(tKsize);
		RHSsize[t].push_back(tRHSsize);

		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			tKsize = 0; tRHSsize = 0;
			if (_mesh.boundaryRegions[r]->dimension) {
				for (auto e = _mesh.boundaryRegions[r]->procNodes->begin(t); e != _mesh.boundaryRegions[r]->procNodes->end(t); ++e) {
					tRHSsize += e->size() * _DOFs;
					tKsize += getMatrixSize(e->size() * _DOFs, mtype);;
				}
			}

			Ksize[t].push_back(tKsize);
			RHSsize[t].push_back(tRHSsize);
		}
	}

	for (size_t i = 0; i < Ksize.front().size(); i++) {
		for (size_t t = 0; t < threads; t++) {
			eslocal tmp = Ksize[t][i];
			Ksize[t][i] = _localKOffset;
			_localKOffset += tmp;

			tmp = RHSsize[t][i];
			RHSsize[t][i] = _localRHSOffset;
			_localRHSOffset += tmp;
		}
	}

	std::vector<IJ> KPattern(_localKOffset);
	std::vector<eslocal> RHSPattern(_localRHSOffset);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		IJ *Koffset = KPattern.data() + Ksize[t][0];
		eslocal *RHSoffset = RHSPattern.data() + RHSsize[t][0];

		auto insert = [&] (serializededata<eslocal, eslocal>::const_iterator &e) {
			eslocal *_RHS = RHSoffset;
			for (auto n = e->begin(); n != e->end(); ++n, ++RHSoffset) {
				*RHSoffset = _DOFMap->datatarray()[*n];
			}
			for (eslocal dof = 1; dof < _DOFs; ++dof) {
				for (size_t n = 0; n < e->size(); ++n, ++RHSoffset) {
					*RHSoffset = *(_RHS + n) + dof;
				}
			}
			insertKPattern(Koffset, _RHS, RHSoffset, mtype);
		};

		for (auto e = _mesh.elements->procNodes->cbegin(t); e != _mesh.elements->procNodes->cend(t); ++e) {
			insert(e);
			Koffset += getMatrixSize(e->size() * _DOFs, mtype);
		}

		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			if (_mesh.boundaryRegions[r]->dimension) {
				for (auto e = _mesh.boundaryRegions[r]->procNodes->cbegin(t); e != _mesh.boundaryRegions[r]->procNodes->cend(t); ++e) {
					insert(e);
					Koffset += getMatrixSize(e->size() * _DOFs, mtype);
				}
			}
		}
	}

	std::vector<eslocal> pK(KPattern.size());
	std::iota(pK.begin(), pK.end(), 0);
	std::sort(pK.begin(), pK.end(), [&] (eslocal i, eslocal j) {
		return KPattern[i] < KPattern[j];
	});

	std::vector<eslocal> pRHS(RHSPattern.size());
	std::iota(pRHS.begin(), pRHS.end(), 0);
	std::sort(pRHS.begin(), pRHS.end(), [&] (eslocal i, eslocal j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	std::vector<std::vector<IJ> > sKBuffer(_mesh.neighbours.size()), rKBuffer(_mesh.neighbours.size());
	std::vector<std::vector<eslocal> > sRHSBuffer(_mesh.neighbours.size()), rRHSBuffer(_mesh.neighbours.size());

	std::vector<eslocal> ndistribution = _mesh.nodes->gatherUniqueNodeDistribution();
	for (size_t n = 0; n < ndistribution.size(); ++n) {
		ndistribution[n] *= _DOFs;
	}
	auto iK = pK.begin();
	auto iRHS = pRHS.begin();
	for (size_t n = 0; n < _mesh.neighbours.size(); ++n) {
		while (KPattern[*iK].row + _mesh.nodes->uniqueOffset < ndistribution[_mesh.neighbours[n]]) {
			if (iK == pK.begin() || KPattern[*iK] != KPattern[*(iK - 1)]) {
				sKBuffer[n].push_back(KPattern[*iK++]);
			}
		}
		while (RHSPattern[*iRHS] + _mesh.nodes->uniqueOffset < ndistribution[_mesh.neighbours[n]]) {
			if (iRHS == pRHS.begin() || RHSPattern[*iRHS] != RHSPattern[*(iRHS - 1)]) {
				sRHSBuffer[n].push_back(RHSPattern[*iRHS++]);
			}
		}
	}

	if (!Communication::receiveUpperUnknownSize(sKBuffer, rKBuffer, _mesh.neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange K pattern.";
	}
	if (!Communication::receiveUpperUnknownSize(sRHSBuffer, rRHSBuffer, _mesh.neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange RHS pattern.";
	}

	for (size_t i = 0; i < rKBuffer.size(); i++) {
		KPattern.insert(KPattern.end(), rKBuffer[i].begin(), rKBuffer[i].end());
	}
	for (size_t i = 0; i < rRHSBuffer.size(); i++) {
		RHSPattern.insert(RHSPattern.end(), rRHSBuffer[i].begin(), rRHSBuffer[i].end());
	}

	size_t localK = pK.size(), localRHS = pRHS.size();
	pK.resize(KPattern.size());
	pRHS.resize(RHSPattern.size());
	std::iota(pK.begin() + localK, pK.end(), localK);
	std::iota(pRHS.begin() + localRHS, pRHS.end(), localRHS);
	std::sort(pK.begin(), pK.end(), [&] (eslocal i, eslocal j) {
		return KPattern[i] < KPattern[j];
	});
	std::sort(pRHS.begin(), pRHS.end(), [&] (eslocal i, eslocal j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	_KPermutation.resize(KPattern.size());
	_RHSPermutation.resize(RHSPattern.size());
	_instance.K.resize(1);
	_instance.M.resize(1);
	_instance.f.resize(1);
	_instance.R.resize(1);
	_instance.primalSolution.resize(1);
	_instance.K.front().haloRows = (_mesh.nodes->size - _mesh.nodes->uniqueSize) * _DOFs;
	_instance.K.front().rows = _mesh.nodes->size * _DOFs;
	_instance.K.front().cols = _mesh.nodes->uniqueTotalSize * _DOFs;
	_instance.f.front().resize(_mesh.nodes->size * _DOFs);
	_instance.R.front().resize(_mesh.nodes->size * _DOFs);
	_instance.primalSolution.front().resize(_mesh.nodes->size * _DOFs);

	std::vector<eslocal> &ROW = _instance.K.front().CSR_I_row_indices;
	std::vector<eslocal> &COL = _instance.K.front().CSR_J_col_indices;
	std::vector<double> &VAL  = _instance.K.front().CSR_V_values;

	ROW.reserve(_mesh.nodes->size * _DOFs + 1);
	ROW.push_back(1);
	COL.push_back(KPattern[pK.front()].column + 1);
	_KPermutation[pK.front()] = 0;
	_RHSPermutation[pRHS.front()] = 0;
	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
			++nonzeros;
			COL.push_back(KPattern[pK[i]].column + 1);
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

	_instance.K.front().nnz = COL.size();
	_instance.K.front().mtype = mtype;
	switch (mtype) {
	case MatrixType::REAL_UNSYMMETRIC: _instance.K.front().type = 'G'; break;
	default: _instance.K.front().type = 'S';
	}
	_instance.K.front().CSR_I_row_indices.swap(ROW);
	_instance.K.front().CSR_J_col_indices.swap(COL);
	_instance.M.front() = _instance.K.front();
	_instance.M.front().mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	_instance.M.front().type = 'S';
}

void UniformNodesComposer::assemble(Matrices matrices)
{
//	MatrixType mtype = _controler.getMatrixType();
	MatrixType mtype = MatrixType::REAL_UNSYMMETRIC; // HYPRE not support symmetric systems

	_controler.updateData();
	clearMatrices(matrices, 0);

	#pragma omp parallel for
	for  (size_t d = 0; d < _instance.domains; d++) {
		size_t KIndex = 0, RHSIndex = 0;
		double KReduction = 1, RHSReduction = _step.internalForceReduction;
		Controler::InstanceFiller filler;

		switch (mtype) {
		case MatrixType::REAL_UNSYMMETRIC:
			filler.insert = [&] (size_t size) {
				for (auto r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						#pragma omp atomic
						_instance.f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						#pragma omp atomic
						_instance.R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (auto c = 0; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							#pragma omp atomic
							_instance.K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							#pragma omp atomic
							_instance.M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						_instance.f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						_instance.R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							_instance.K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							_instance.M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		}

		filler.begin = _mesh.elements->elementsDistribution[d];
		filler.end = _mesh.elements->elementsDistribution[d + 1];

		_controler.processElements(matrices, filler);

//		KReduction = _step.internalForceReduction;
//
//		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
//			if (_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {
//				filler.begin = _mesh.elements->elementsDistribution[d];
//				filler.end = _mesh.elements->elementsDistribution[d + 1];
//			}
//		}

//		auto boundary = [&] (Matrices restriction) {
//			for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
//				if (_mesh.boundaryRegions[r]->dimension == 2) {
//					if (_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {
//						eslocal begin = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d]].begin;
//						eslocal end = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
//						auto nodes = _mesh.boundaryRegions[r]->procNodes->cbegin() + begin;
//						for (eslocal i = begin; i < end; ++i, ++nodes) {
//							processFace(d, _mesh.boundaryRegions[r], restriction, i, Ke, Me, Re, fe);
//							switch (getMatrixType(domain)) {
//							case MatrixType::REAL_UNSYMMETRIC: fullInsert(nodes, _step.internalForceReduction); break;
//							default: upperInsert(nodes, _step.internalForceReduction);
//							}
//						}
//					}
//				}
//				if (_mesh.boundaryRegions[r]->dimension == 1) {
//					if (_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {
//						eslocal begin = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d]].begin;
//						eslocal end = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
//						auto nodes = _mesh.boundaryRegions[r]->procNodes->cbegin() + begin;
//						for (eslocal i = begin; i < end; ++i, ++nodes) {
//							processEdge(d, _mesh.boundaryRegions[r], restriction, i, Ke, Me, Re, fe);
//							switch (getMatrixType(domain)) {
//							case MatrixType::REAL_UNSYMMETRIC: fullInsert(nodes, _step.internalForceReduction); break;
//							default: upperInsert(nodes, _step.internalForceReduction);
//							}
//						}
//					}
//				}
//			}
//		};
	}
}

void UniformNodesComposer::setDirichlet()
{
	auto &ROW = _instance.K.front().CSR_I_row_indices;
	auto &COL = _instance.K.front().CSR_J_col_indices;
	auto &VAL = _instance.K.front().CSR_V_values;
	auto &RHS = _instance.f.front();

	std::vector<double> values(_dirichletMap.size());
	_controler.dirichletValues(values);

	for (size_t i = 0; i < _dirichletMap.size(); ++i) {
		RHS[_dirichletMap[i]] = values[_dirichletPermutation[i]];
		for (eslocal j = ROW[_dirichletMap[i]]; j < ROW[_dirichletMap[i]+ 1]; j++) {
			if (COL[j - 1] - 1 == _dirichletMap[i]) {
				VAL[j - 1] = 1;
			} else {
				VAL[j - 1] = 0;
				eslocal r = COL[j - 1] - 1 - _mesh.nodes->uniqueOffset;
				for (eslocal c = ROW[r]; c < ROW[r + 1]; c++) {
					if (COL[c - 1] - 1 == _dirichletMap[i] + _mesh.nodes->uniqueOffset) {
						RHS[r] -= VAL[c - 1] * RHS[_dirichletMap[i]];
						VAL[c - 1] = 0;
					}
				}
			}
		}
	}
}

void UniformNodesComposer::synchronize()
{
	std::vector<std::vector<double> > sBuffer(_mesh.neighbours.size()), rBuffer(_mesh.neighbours.size());

	size_t n = 0;
	for (auto it = _mesh.nodes->pintervals.begin(); it != _mesh.nodes->pintervals.end() && it->sourceProcess < environment->MPIrank; ++it) {
		if (_mesh.neighbours[n] < it->sourceProcess) {
			++n;
		}
		sBuffer[n].push_back(it->end - it->begin);
		sBuffer[n].insert(sBuffer[n].end(), _instance.f.front().begin() + it->begin, _instance.f.front().begin() + it->end);
	}

	n = 0;
	for (auto it = _mesh.nodes->pintervals.begin(); it != _mesh.nodes->pintervals.end() && it->sourceProcess < environment->MPIrank; ++it) {
		if (_mesh.neighbours[n] < it->sourceProcess) {
			++n;
		}
		auto begin = _instance.K.front().CSR_I_row_indices[it->begin] - 1;
		auto end = _instance.K.front().CSR_I_row_indices[it->end] - 1;
		sBuffer[n].insert(sBuffer[n].end(), _instance.K.front().CSR_V_values.begin() + begin, _instance.K.front().CSR_V_values.begin() + end);
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, _mesh.neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange CSR pattern.";
	}

	size_t KIndex = _localKOffset, RHSIndex = _localRHSOffset;
	for (size_t i = 0, j = 1; i < rBuffer.size(); ++i, j = 0) {
		for (; j < rBuffer[i].front(); ++j, ++RHSIndex) {
			_instance.f.front()[_RHSPermutation[RHSIndex]] += rBuffer[i][j];
		}
		for (; j < rBuffer[i].size(); ++j, ++KIndex) {
			_instance.K.front().CSR_V_values[_KPermutation[KIndex]] += rBuffer[i][j];
		}
	}
}
