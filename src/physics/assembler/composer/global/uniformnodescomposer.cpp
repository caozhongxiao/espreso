

#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "uniformnodescomposer.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
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
	MatrixType mtype = _provider.getMatrixType();

	_nDistribution = info::mesh->nodes->gatherUniqueNodeDistribution();
	for (size_t n = 0; n < _nDistribution.size(); ++n) {
		_nDistribution[n] *= _DOFs;
	}

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

//	for (size_t i = 0; i < KPattern.size(); i++) {
//		printf("[%d,%d] ", KPattern[i].row, KPattern[i].column);
//	}
//	printf("\n");

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
	std::sort(pK.begin(), pK.end(), [&] (esint i, esint j) {
		return KPattern[i] < KPattern[j];
	});
	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
		return RHSPattern[i] < RHSPattern[j];
	});

//	for (size_t i = 0; i < KPattern.size(); i++) {
//		std::cout << KPattern[i].row << ":" << KPattern[i].column << "\n";
//	}
//	std::cout << RHSPattern;

	_KPermutation.resize(KPattern.size());
	_RHSPermutation.resize(RHSPattern.size());
	data->K.resize(1);
	data->origK.resize(1);
	data->M.resize(1);
	data->f.resize(1);
	data->R.resize(1);
	data->primalSolution.resize(1);
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

	data->K.front().nnz = COL.size();
	data->K.front().mtype = mtype;
	switch (mtype) {
	case MatrixType::REAL_UNSYMMETRIC: data->K.front().type = 'G'; break;
	default: data->K.front().type = 'S';
	}
	data->K.front().CSR_I_row_indices.swap(ROW);
	data->K.front().CSR_J_col_indices.swap(COL);
	data->M.front() = data->K.front();
	data->M.front().mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	data->M.front().type = 'S';
}

void UniformNodesComposer::assemble(Matrices matrices, const SolverParameters &parameters)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	MatrixType mtype = _provider.getMatrixType();

	clearMatrices(matrices, 0);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t KIndex = _tKOffsets[t], RHSIndex = _tRHSOffsets[t];
		double KReduction = parameters.timeIntegrationConstantK, RHSReduction = parameters.internalForceReduction;
		Controler::InstanceFiller filler;

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
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							#pragma omp atomic
							data->M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						data->f[0][_RHSPermutation[RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						data->R[0][_RHSPermutation[RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							data->K[0].CSR_V_values[_KPermutation[KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							data->M[0].CSR_V_values[_KPermutation[KIndex]] += filler.Me(r, c);
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
}

void UniformNodesComposer::setDirichlet()
{
	auto &ROW = data->K.front().CSR_I_row_indices;
	auto &COL = data->K.front().CSR_J_col_indices;
	auto &VAL = data->K.front().CSR_V_values;
	auto &RHS = data->f.front();

	std::vector<double> values(_dirichletMap.size());
	_controler.dirichletValues(values);

//	std::vector<esint> RROW;
//	for (size_t r = 0; r < data->K.front().CSR_I_row_indices.size() - 1; r++) {
//		RROW.insert(RROW.end(), data->K.front().CSR_I_row_indices[r + 1] - data->K.front().CSR_I_row_indices[r], _DOFMap->datatarray()[r] + 1);
//	}
//
//	int rank = 3;
//
//	Communication::serialize([&] () {
//		if (info::mpi::MPIrank != rank) {
//			return;
//		}
//		printf(" // %d \\\\ \n", info::mpi::MPIrank);
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
//		for (esint r = 0, i = 0, f = 0; r < info::mesh->nodes->uniqueTotalSize; r++) {
//			for (esint c = 0; c < info::mesh->nodes->uniqueTotalSize; c++) {
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
//		std::cout << data->K.front();
//		std::cout << data->f.front();
//
//		for (size_t i = 0; i < _dirichletMap.size(); i++) {
//			std::cout << _dirichletMap[i] + 1 << " ";
//		}
//		std::cout << "\n";
//	});

	for (size_t i = 0; i < _dirichletMap.size(); ++i) {
		RHS[_dirichletMap[i]] = values[_dirichletPermutation[i]];
//		if (info::mpi::MPIrank == rank) {
//			std::cout << "RHS[" << _DOFMap->datatarray()[_dirichletMap[i]] + 1 << "] = " << values[_dirichletPermutation[i]] << "\n";
//		}
		esint col = _DOFMap->datatarray()[_dirichletMap[i]] + 1;
		for (esint j = ROW[_dirichletMap[i]]; j < ROW[_dirichletMap[i] + 1]; j++) {
			if (COL[j - 1] == col) {
//				if (info::mpi::MPIrank == rank) {
//					std::cout << "[" << _DOFMap->datatarray()[_dirichletMap[i]] + 1 << ":" << COL[j - 1] << "] = 1\n";
//				}
				VAL[j - 1] = 1;
			} else {
//				if (info::mpi::MPIrank == rank) {
//					std::cout << "[" << _DOFMap->datatarray()[_dirichletMap[i]] + 1 << ":" << COL[j - 1] << "] = 0\n";
//				}
				VAL[j - 1] = 0;
				size_t r = std::lower_bound(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end(), COL[j - 1] - 1) - _DOFMap->datatarray().begin();
				if (r < _DOFMap->datatarray().size() && _DOFMap->datatarray()[r] == COL[j - 1] - 1) {
//					if (info::mpi::MPIrank == rank) {
//						std::cout << COL[j - 1] << " into " << r << "\n";
//					}
					for (esint c = ROW[r]; c < ROW[r + 1]; c++) {
						if (COL[c - 1] == col) {
//							if (info::mpi::MPIrank == rank) {
//								std::cout << "[" << _DOFMap->datatarray()[r] + 1 << ":" << COL[c - 1] << "] = 0; RHS[" << _DOFMap->datatarray()[r] + 1 << "] -= " << VAL[c - 1] << " * " << RHS[_dirichletMap[i]] << "\n";
//							}
							RHS[r] -= VAL[c - 1] * RHS[_dirichletMap[i]];
							VAL[c - 1] = 0;
						}
					}
				}
			}

//			if (info::mpi::MPIrank == rank) {
//				for (esint r = 0, i = 0, f = 0; r < info::mesh->nodes->uniqueTotalSize; r++) {
//					for (esint c = 0; c < info::mesh->nodes->uniqueTotalSize; c++) {
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
//		if (info::mpi::MPIrank != rank) {
//			return;
//		}
//		for (size_t i = 0; i < _dirichletMap.size(); i++) {
//			printf("%d ", _dirichletMap[i]);
//		}
//		printf("\n");
//		printf(" // %d \\\\ \n", info::mpi::MPIrank);
//		for (esint r = 0, i = 0, f = 0; r < info::mesh->nodes->uniqueTotalSize; r++) {
//			for (esint c = 0; c < info::mesh->nodes->uniqueTotalSize; c++) {
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
//		std::cout << data->K.front();
//	});
}

void UniformNodesComposer::synchronize()
{
	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	auto nranks = info::mesh->nodes->ranks->begin();
	auto DOFs = _DOFMap->begin();
	for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
		esint r = 0;
		while (info::mesh->neighbours[r] < nranks->front()) {
			++r;
		}
		for (int dof = 0; dof < _DOFs; ++dof) {
			sBuffer[r].push_back(data->f.front()[n * _DOFs + dof]);
		}
		_nRHSSize[r] = sBuffer[r].size();
	}

	nranks = info::mesh->nodes->ranks->begin();
	DOFs = _DOFMap->begin();
	for (esint n = 0; n < info::mesh->nodes->size && DOFs->front() < _nDistribution[info::mpi::rank]; ++n, ++nranks, ++DOFs) {
		esint r = 0;
		while (info::mesh->neighbours[r] < nranks->front()) {
			++r;
		}
		auto begin = data->K.front().CSR_I_row_indices[DOFs->begin() - _DOFMap->datatarray().begin()] - 1;
		auto end = data->K.front().CSR_I_row_indices[DOFs->end() - _DOFMap->datatarray().begin()] - 1;
		sBuffer[r].insert(sBuffer[r].end(), data->K.front().CSR_V_values.begin() + begin, data->K.front().CSR_V_values.begin() + end);
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize assembled data.";
	}

	size_t KIndex = _localKOffset, RHSIndex = _localRHSOffset;
	for (size_t i = 0, j = 0; i < rBuffer.size(); ++i, j = 0) {
		for (; rBuffer[i].size() && j < (size_t)_nRHSSize[i]; ++j, ++RHSIndex) {
			data->f.front()[_RHSPermutation[RHSIndex]] += rBuffer[i][j];
		}
		for (; j < rBuffer[i].size(); ++j, ++KIndex) {
			data->K.front().CSR_V_values[_KPermutation[KIndex]] += rBuffer[i][j];
		}
	}

//	std::cout << data->K.front();
//	std::cout << data->f.front();
}

void UniformNodesComposer::fillSolution()
{
	std::vector<double> &solution = _controler.solution()->data;

	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	size_t RHSIndex = _localRHSOffset;
	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		if (info::mesh->neighbours[n] < info::mpi::rank) {
			rBuffer[n].resize(_nRHSSize[n]);
		} else {
			sBuffer[n].reserve(_nRHSSize[n]);
			for (esint i = 0; i < _nRHSSize[n]; ++i) {
				sBuffer[n].push_back(data->primalSolution.front()[_RHSPermutation[RHSIndex++]]);
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
			solution[n * _DOFs + dof] = rBuffer[r][rIndices[r]++];
		}
	}

	memcpy(
			solution.data() + data->K.front().haloRows,
			data->primalSolution.front().data() + data->K.front().haloRows,
			sizeof(double) * (data->K.front().rows - data->K.front().haloRows));
}
