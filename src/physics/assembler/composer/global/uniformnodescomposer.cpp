
#include "uniformnodescomposer.h"

#include "../../controllers/controller.h"

#include "../../../../globals/run.h"
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
#include "../../../dataholder.h"

using namespace espreso;

void UniformNodesComposer::initDOFs()
{
	// ASSUME THAT SHARED NODES ARE SORTED IN THE SAME ORDER ON ALL PROCESSES

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> doffset(threads);
	std::vector<std::vector<eslocal> > roffset(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = run::mesh->nodes->ranks->begin(t);
		eslocal dsize = 0;
		std::vector<eslocal> troffset(run::mesh->neighbours.size());

		for (size_t n = run::mesh->nodes->distribution[t]; n < run::mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				dsize += _DOFs;
			} else {
				eslocal noffset = 0;
				while (run::mesh->neighbours[noffset] < ranks->front()) {
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

	std::vector<std::vector<std::vector<eslocal> > > sBuffer(threads, std::vector<std::vector<eslocal> >(run::mesh->neighbours.size()));
	std::vector<std::vector<eslocal> > rBuffer(run::mesh->neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = run::mesh->nodes->ranks->begin(t);
		eslocal toffset = doffset[t];
		std::vector<std::vector<eslocal> > tBuffer(run::mesh->neighbours.size());

		for (size_t n = run::mesh->nodes->distribution[t]; n < run::mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				eslocal noffset = 0;
				for (auto r = ranks->begin() + 1; r != ranks->end(); ++r) {
					while (run::mesh->neighbours[noffset] < *r) {
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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	std::vector<std::vector<eslocal> > DOFDistribution(threads), DOFData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = run::mesh->nodes->ranks->begin(t);
		eslocal toffset = doffset[t];
		std::vector<eslocal> troffset = roffset[t];
		std::vector<eslocal> tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}

		for (size_t n = run::mesh->nodes->distribution[t]; n < run::mesh->nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				for (int dof = 0; dof < _DOFs; ++dof) {
					tdata.push_back(goffset + toffset++);
				}
			} else {
				eslocal noffset = 0;
				while (run::mesh->neighbours[noffset] < ranks->front()) {
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

	_nDistribution = run::mesh->nodes->gatherUniqueNodeDistribution();
	for (size_t n = 0; n < _nDistribution.size(); ++n) {
		_nDistribution[n] *= _DOFs;
	}

	_tKOffsets.resize(threads);
	_tRHSOffsets.resize(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		eslocal tKsize = 0, tRHSsize = 0;
		for (auto e = run::mesh->elements->procNodes->begin(t); e != run::mesh->elements->procNodes->end(t); ++e) {
			tRHSsize += e->size() * _DOFs;
			tKsize += getMatrixSize(e->size() * _DOFs, mtype);
		}

		for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
			if (run::mesh->boundaryRegions[r]->dimension) {
				for (auto e = run::mesh->boundaryRegions[r]->procNodes->begin(t); e != run::mesh->boundaryRegions[r]->procNodes->end(t); ++e) {
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
	std::vector<eslocal> RHSPattern(_localRHSOffset);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		IJ *Koffset = KPattern.data() + _tKOffsets[t];
		eslocal *RHSoffset = RHSPattern.data() + _tRHSOffsets[t];

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

		for (auto e = run::mesh->elements->procNodes->cbegin(t); e != run::mesh->elements->procNodes->cend(t); ++e) {
			insert(e);
			Koffset += getMatrixSize(e->size() * _DOFs, mtype);
		}

		for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
			if (run::mesh->boundaryRegions[r]->dimension) {
				for (auto e = run::mesh->boundaryRegions[r]->procNodes->cbegin(t); e != run::mesh->boundaryRegions[r]->procNodes->cend(t); ++e) {
					insert(e);
					Koffset += getMatrixSize(e->size() * _DOFs, mtype);
				}
			}
		}
	}

//	for (size_t i = 0; i < KPattern.size(); i++) {
//		printf("[%d,%d] ", KPattern[i].row, KPattern[i].column);
//	}
//	printf("\n");

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

	std::vector<std::vector<IJ> > sKBuffer(run::mesh->neighbours.size()), rKBuffer(run::mesh->neighbours.size());
	std::vector<std::vector<eslocal> > sRHSBuffer(run::mesh->neighbours.size()), rRHSBuffer(run::mesh->neighbours.size());

	auto iK = pK.begin();
	auto iRHS = pRHS.begin();
	for (size_t n = 0; n < run::mesh->neighbours.size() && run::mesh->neighbours[n] < environment->MPIrank; ++n) {
		while (KPattern[*iK].row < _nDistribution[run::mesh->neighbours[n] + 1]) {
			if (iK == pK.begin() || KPattern[*iK] != KPattern[*(iK - 1)]) {
				sKBuffer[n].push_back(KPattern[*iK]);
			}
			++iK;
		}
		while (RHSPattern[*iRHS] < _nDistribution[run::mesh->neighbours[n] + 1]) {
			if (iRHS == pRHS.begin() || RHSPattern[*iRHS] != RHSPattern[*(iRHS - 1)]) {
				sRHSBuffer[n].push_back(RHSPattern[*iRHS]);
			}
			++iRHS;
		}
	}

	if (!Communication::receiveUpperUnknownSize(sKBuffer, rKBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange K pattern.";
	}
	if (!Communication::receiveUpperUnknownSize(sRHSBuffer, rRHSBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange RHS pattern.";
	}

	_nKSize.clear();
	_nRHSSize.clear();
	for (size_t i = 0; i < rKBuffer.size(); i++) {
		KPattern.insert(KPattern.end(), rKBuffer[i].begin(), rKBuffer[i].end());
		_nKSize.push_back(rKBuffer[i].size());
	}
	for (size_t i = 0; i < rRHSBuffer.size(); i++) {
		RHSPattern.insert(RHSPattern.end(), rRHSBuffer[i].begin(), rRHSBuffer[i].end());
		_nRHSSize.push_back(rRHSBuffer[i].size());
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

//	std::cout << KPattern;
//	std::cout << RHSPattern;

	_KPermutation.resize(KPattern.size());
	_RHSPermutation.resize(RHSPattern.size());
	run::data->K.resize(1);
	run::data->M.resize(1);
	run::data->f.resize(1);
	run::data->R.resize(1);
	run::data->primalSolution.resize(1);
	run::data->K.front().haloRows = (run::mesh->nodes->size - run::mesh->nodes->uniqueSize) * _DOFs;
	run::data->K.front().rows = run::mesh->nodes->size * _DOFs;
	run::data->K.front().cols = run::mesh->nodes->uniqueTotalSize * _DOFs;
	run::data->f.front().resize(run::mesh->nodes->size * _DOFs);
	run::data->R.front().resize(run::mesh->nodes->size * _DOFs);
	run::data->primalSolution.front().resize(run::mesh->nodes->size * _DOFs);

	std::vector<eslocal> &ROW = run::data->K.front().CSR_I_row_indices;
	std::vector<eslocal> &COL = run::data->K.front().CSR_J_col_indices;
	std::vector<double> &VAL  = run::data->K.front().CSR_V_values;

	ROW.reserve(run::mesh->nodes->size * _DOFs + 1);
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

	run::data->K.front().nnz = COL.size();
	run::data->K.front().mtype = mtype;
	switch (mtype) {
	case MatrixType::REAL_UNSYMMETRIC: run::data->K.front().type = 'G'; break;
	default: run::data->K.front().type = 'S';
	}
	run::data->K.front().CSR_I_row_indices.swap(ROW);
	run::data->K.front().CSR_J_col_indices.swap(COL);
	run::data->M.front() = run::data->K.front();
	run::data->M.front().mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	run::data->M.front().type = 'S';
}

void UniformNodesComposer::assemble(Matrices matrices)
{
	_controler.nextTime();

	size_t threads = environment->OMP_NUM_THREADS;

//	MatrixType mtype = _controler.getMatrixType();
	MatrixType mtype = MatrixType::REAL_UNSYMMETRIC; // HYPRE not support symmetric systems

	clearMatrices(matrices, 0);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t KIndex = _tKOffsets[t], RHSIndex = _tRHSOffsets[t];
		double KReduction = 1, RHSReduction = 1; //_step.internalForceReduction;
		Controler::InstanceFiller filler;

		switch (mtype) {
		case MatrixType::REAL_UNSYMMETRIC:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						#pragma omp atomic
						run::data->f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						#pragma omp atomic
						run::data->R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = 0; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							#pragma omp atomic
							run::data->K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							#pragma omp atomic
							run::data->M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						run::data->f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						run::data->R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							run::data->K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							run::data->M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		}

		filler.begin = run::mesh->elements->distribution[t];
		filler.end = run::mesh->elements->distribution[t + 1];

		_controler.processElements(matrices, filler);

		KReduction = 1; // _step.internalForceReduction;

		for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
			if (run::mesh->boundaryRegions[r]->distribution.size()) {
				filler.begin = run::mesh->boundaryRegions[r]->distribution[t];
				filler.end = run::mesh->boundaryRegions[r]->distribution[t + 1];
				_controler.processBoundary(matrices, r, filler);
			}
		}
	}
}

void UniformNodesComposer::setDirichlet()
{
	auto &ROW = run::data->K.front().CSR_I_row_indices;
	auto &COL = run::data->K.front().CSR_J_col_indices;
	auto &VAL = run::data->K.front().CSR_V_values;
	auto &RHS = run::data->f.front();

	std::vector<double> values(_dirichletMap.size());
	_controler.dirichletValues(values);

//	std::vector<eslocal> RROW;
//	for (size_t r = 0; r < run::data->K.front().CSR_I_row_indices.size() - 1; r++) {
//		RROW.insert(RROW.end(), run::data->K.front().CSR_I_row_indices[r + 1] - run::data->K.front().CSR_I_row_indices[r], _DOFMap->datatarray()[r] + 1);
//	}
//
//	int rank = 3;
//
//	Communication::serialize([&] () {
//		if (environment->MPIrank != rank) {
//			return;
//		}
//		printf(" // %d \\\\ \n", environment->MPIrank);
//		for (size_t i = 0; i < _dirichletMap.size(); i++) {
//			printf("%d ", _DOFMap->datatarray()[_dirichletMap[i]] + 1);
//		}
//		printf("\n");
//
//		for (size_t i = 0; i < _DOFMap->datatarray().size(); i++) {
//			printf("%d ", _DOFMap->datatarray()[i] + 1);
//		}
//		printf("\n");
//
//		for (eslocal r = 0, i = 0, f = 0; r < run::mesh->nodes->uniqueTotalSize; r++) {
//			for (eslocal c = 0; c < run::mesh->nodes->uniqueTotalSize; c++) {
//				if (i < RROW.size() && RROW[i] == r + 1 && COL[i] == c + 1) {
//					if (VAL[i] > -0.00001) {
//						if (std::fabs(VAL[i]) > 10) {
//							printf(" %3.1f ", VAL[i++]);
//						} else {
//							printf(" %3.2f ", VAL[i++]);
//						}
//					} else {
//						if (std::fabs(VAL[i]) > 10) {
//							printf("%3.1f ", VAL[i++]);
//						} else {
//							printf("%3.2f ", VAL[i++]);
//						}
//					}
//				} else {
//					printf("      ");
//				}
//			}
//			if (f < _DOFMap->datatarray().size() && _DOFMap->datatarray()[f] == r) {
//				printf(" = %3.2f\n", RHS[f++]);
//			} else {
//				printf(" =\n");
//			}
//		}
//		printf("------------------\n");
//	});
//
//	Communication::serialize([&] () {
//		std::cout << run::data->K.front();
//		std::cout << run::data->f.front();
//
//		for (size_t i = 0; i < _dirichletMap.size(); i++) {
//			std::cout << _dirichletMap[i] + 1 << " ";
//		}
//		std::cout << "\n";
//	});

	auto ndofbegin = _DOFMap->datatarray().begin();
	auto ndofend = (_DOFMap->begin() + (run::mesh->nodes->size - run::mesh->nodes->uniqueSize))->begin();

	for (size_t i = 0; i < _dirichletMap.size(); ++i) {
		RHS[_dirichletMap[i]] = values[_dirichletPermutation[i]];
//		if (environment->MPIrank == rank) {
//			std::cout << "RHS[" << _DOFMap->datatarray()[_dirichletMap[i]] + 1 << "] = " << values[_dirichletPermutation[i]] << "\n";
//		}
		eslocal col = _DOFMap->datatarray()[_dirichletMap[i]] + 1;
		for (eslocal j = ROW[_dirichletMap[i]]; j < ROW[_dirichletMap[i] + 1]; j++) {
			if (COL[j - 1] == col) {
//				if (environment->MPIrank == rank) {
//					std::cout << "[" << _DOFMap->datatarray()[_dirichletMap[i]] + 1 << ":" << COL[j - 1] << "] = 1\n";
//				}
				VAL[j - 1] = 1;
			} else {
//				if (environment->MPIrank == rank) {
//					std::cout << "[" << _DOFMap->datatarray()[_dirichletMap[i]] + 1 << ":" << COL[j - 1] << "] = 0\n";
//				}
				VAL[j - 1] = 0;
				eslocal r = std::lower_bound(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end(), COL[j - 1] - 1) - _DOFMap->datatarray().begin();
				if (r < _DOFMap->datatarray().size() && _DOFMap->datatarray()[r] == COL[j - 1] - 1) {
//					if (environment->MPIrank == rank) {
//						std::cout << COL[j - 1] << " into " << r << "\n";
//					}
					for (eslocal c = ROW[r]; c < ROW[r + 1]; c++) {
						if (COL[c - 1] == col) {
//							if (environment->MPIrank == rank) {
//								std::cout << "[" << _DOFMap->datatarray()[r] + 1 << ":" << COL[c - 1] << "] = 0; RHS[" << _DOFMap->datatarray()[r] + 1 << "] -= " << VAL[c - 1] << " * " << RHS[_dirichletMap[i]] << "\n";
//							}
							RHS[r] -= VAL[c - 1] * RHS[_dirichletMap[i]];
							VAL[c - 1] = 0;
						}
					}
				}
			}

//			if (environment->MPIrank == rank) {
//				for (eslocal r = 0, i = 0, f = 0; r < run::mesh->nodes->uniqueTotalSize; r++) {
//					for (eslocal c = 0; c < run::mesh->nodes->uniqueTotalSize; c++) {
//						if (i < RROW.size() && RROW[i] == r + 1 && COL[i] == c + 1) {
//							if (VAL[i] > -0.00001) {
//								if (std::fabs(VAL[i]) > 10) {
//									printf(" %3.1f ", VAL[i++]);
//								} else {
//									printf(" %3.2f ", VAL[i++]);
//								}
//							} else {
//								if (std::fabs(VAL[i]) > 10) {
//									printf("%3.1f ", VAL[i++]);
//								} else {
//									printf("%3.2f ", VAL[i++]);
//								}
//							}
//						} else {
//							printf("      ");
//						}
//					}
//					if (f < _DOFMap->datatarray().size() && _DOFMap->datatarray()[f] == r) {
//						printf(" = %3.2f\n", RHS[f++]);
//					} else {
//						printf(" =\n");
//					}
//				}
//				printf("------------------\n");
//			}
		}
	}


//	Communication::serialize([&] () {
//		if (environment->MPIrank != rank) {
//			return;
//		}
//		for (size_t i = 0; i < _dirichletMap.size(); i++) {
//			printf("%d ", _dirichletMap[i]);
//		}
//		printf("\n");
//		printf(" // %d \\\\ \n", environment->MPIrank);
//		for (eslocal r = 0, i = 0, f = 0; r < run::mesh->nodes->uniqueTotalSize; r++) {
//			for (eslocal c = 0; c < run::mesh->nodes->uniqueTotalSize; c++) {
//				if (i < RROW.size() && RROW[i] == r + 1 && COL[i] == c + 1) {
//					if (VAL[i] > -0.00001) {
//						if (std::fabs(VAL[i]) > 10) {
//							printf(" %3.1f ", VAL[i++]);
//						} else {
//							printf(" %3.2f ", VAL[i++]);
//						}
//					} else {
//						if (std::fabs(VAL[i]) > 10) {
//							printf("%3.1f ", VAL[i++]);
//						} else {
//							printf("%3.2f ", VAL[i++]);
//						}
//					}
//				} else {
//					printf("      ");
//				}
//			}
//			if (f < _DOFMap->datatarray().size() && _DOFMap->datatarray()[f] == r) {
//				printf(" = %3.2f\n", RHS[f++]);
//			} else {
//				printf(" =\n");
//			}
//		}
//		printf("------------------\n");
//	});
//	Communication::serialize([&] () {
//		std::cout << run::data->K.front();
//	});
}

void UniformNodesComposer::synchronize()
{
	std::vector<std::vector<double> > sBuffer(run::mesh->neighbours.size()), rBuffer(run::mesh->neighbours.size());

	auto nranks = run::mesh->nodes->ranks->begin();
	auto DOFs = _DOFMap->begin();
	for (eslocal n = 0; n < run::mesh->nodes->size && DOFs->front() < _nDistribution[environment->MPIrank]; ++n, ++nranks, ++DOFs) {
		eslocal r = 0;
		while (run::mesh->neighbours[r] < nranks->front()) {
			++r;
		}
		for (int dof = 0; dof < _DOFs; ++dof) {
			sBuffer[r].push_back(run::data->f.front()[n * _DOFs + dof]);
		}
		_nRHSSize[r] = sBuffer[r].size();
	}

	nranks = run::mesh->nodes->ranks->begin();
	DOFs = _DOFMap->begin();
	for (eslocal n = 0; n < run::mesh->nodes->size && DOFs->front() < _nDistribution[environment->MPIrank]; ++n, ++nranks, ++DOFs) {
		eslocal r = 0;
		while (run::mesh->neighbours[r] < nranks->front()) {
			++r;
		}
		auto begin = run::data->K.front().CSR_I_row_indices[DOFs->begin() - _DOFMap->datatarray().begin()] - 1;
		auto end = run::data->K.front().CSR_I_row_indices[DOFs->end() - _DOFMap->datatarray().begin()] - 1;
		sBuffer[r].insert(sBuffer[r].end(), run::data->K.front().CSR_V_values.begin() + begin, run::data->K.front().CSR_V_values.begin() + end);
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize assembled data.";
	}

	size_t KIndex = _localKOffset, RHSIndex = _localRHSOffset;
	for (size_t i = 0, j = 0; i < rBuffer.size(); ++i, j = 0) {
		for (; rBuffer[i].size() && j < (size_t)_nRHSSize[i]; ++j, ++RHSIndex) {
			run::data->f.front()[_RHSPermutation[RHSIndex]] += rBuffer[i][j];
		}
		for (; j < rBuffer[i].size(); ++j, ++KIndex) {
			run::data->K.front().CSR_V_values[_KPermutation[KIndex]] += rBuffer[i][j];
		}
	}

//	std::cout << run::data->K.front();
//	std::cout << run::data->f.front();
}

void UniformNodesComposer::fillSolution()
{
	std::vector<double> &solution = _controler.getSolutionStore();

	std::vector<std::vector<double> > sBuffer(run::mesh->neighbours.size()), rBuffer(run::mesh->neighbours.size());

	size_t RHSIndex = _localRHSOffset;
	for (size_t n = 0; n < run::mesh->neighbours.size(); ++n) {
		if (run::mesh->neighbours[n] < environment->MPIrank) {
			rBuffer[n].resize(_nRHSSize[n]);
		} else {
			sBuffer[n].reserve(_nRHSSize[n]);
			for (eslocal i = 0; i < _nRHSSize[n]; ++i) {
				sBuffer[n].push_back(run::data->primalSolution.front()[_RHSPermutation[RHSIndex++]]);
			}
		}
	}

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize assembled data.";
	}

	std::vector<eslocal> rIndices(run::mesh->neighbours.size());
	auto nranks = run::mesh->nodes->ranks->begin();
	auto DOFs = _DOFMap->begin();
	for (eslocal n = 0; n < run::mesh->nodes->size && DOFs->front() < _nDistribution[environment->MPIrank]; ++n, ++nranks, ++DOFs) {
		eslocal r = 0;
		while (run::mesh->neighbours[r] < nranks->front()) {
			++r;
		}
		for (int dof = 0; dof < _DOFs; ++dof) {
			solution[n * _DOFs + dof] = rBuffer[r][rIndices[r]++];
		}
	}

	memcpy(
			solution.data() + run::data->K.front().haloRows,
			run::data->primalSolution.front().data() + run::data->K.front().haloRows,
			sizeof(double) * (run::data->K.front().rows - run::data->K.front().haloRows));

//	Communication::serialize([&] () {
//		std::cout << environment->MPIrank << "\n";
//		std::cout << run::data->primalSolution.front();
//		std::cout << solution;
//		std::cout << " >> <<\n";
//	});
}
