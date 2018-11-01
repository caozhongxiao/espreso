
#include "physicsinvectors.h"

#include "../instance.h"

#include "../../basis/utilities/communication.h"
#include "../../mesh/mesh.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/utilities/utils.h"
#include "../../solver/generic/SparseMatrix.h"

#include <numeric>
#include <algorithm>

using namespace espreso;


struct IJ {
	eslocal row, column;
};

bool operator==(const IJ &left, const IJ &right)
{
	return left.row == right.row && left.column == right.column;
}

bool operator!=(const IJ &left, const IJ &right)
{
	return !(left == right);
}

PhysicsInVectors::PhysicsInVectors(const std::string &name, Mesh &mesh, Instance &instance, Step &step, const PhysicsConfiguration &configuration)
: _name(name), _mesh(mesh), _instance(instance), _step(step), _configuration(configuration), _DOFs{0, 0, 0, 1, 1}, _domainNodes(NULL),
  _localK(0), _localRHS(0)
{
	size_t threads = environment->OMP_NUM_THREADS;
	_nDistribution.resize(threads + 1);

	for (size_t t = 0; t < threads; t++) {
		_nDistribution[t + 1] = _nDistribution[t];
		for (eslocal d = _mesh.elements->domainDistribution[t]; d != _mesh.elements->domainDistribution[t + 1]; ++d) {
			eslocal dbegin = _mesh.elements->elementsDistribution[d];
			eslocal dend = _mesh.elements->elementsDistribution[d + 1];
			_nDistribution[t + 1] += (_mesh.elements->procNodes->begin() + dend)->begin() - (_mesh.elements->procNodes->begin() + dbegin)->begin();
		}
	}
}

void PhysicsInVectors::buildCSRPattern()
{
	size_t threads = environment->OMP_NUM_THREADS;

//	for (size_t d = 0; d < _instance.domains; d++) {
//		std::cout << " -- DOMAIN " << d << " -- \n";
//		std::cout << "R: " << _instance.K.front().CSR_I_row_indices;
//		std::cout << "C: " << _instance.K.front().CSR_J_col_indices;
//	}

	_domainNodes = new serializededata<eslocal, eslocal>(*_mesh.elements->procNodes);
	_DOFsPermutation.resize(_mesh.elements->ndomains);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh.elements->domainDistribution[t]; d != _mesh.elements->domainDistribution[t + 1]; ++d) {
			eslocal dbegin = _mesh.elements->elementsDistribution[d];
			eslocal dend = _mesh.elements->elementsDistribution[d + 1];
			std::vector<eslocal> nodes, DOFs;
			std::vector<IJ> pattern;
			nodes.insert(nodes.end(), (_mesh.elements->procNodes->begin() + dbegin)->begin(), (_mesh.elements->procNodes->begin() + dend)->begin());
			Esutils::sortAndRemoveDuplicity(nodes);

//			std::cout << "NODES: " << nodes;

			auto ebegin = _domainNodes->begin() + dbegin;
			auto eend = _domainNodes->begin() + dend;
			for (auto n = ebegin->begin(); n != eend->begin(); ++n) {
				*n = std::lower_bound(nodes.begin(), nodes.end(), *n) - nodes.begin();
			}

			for (auto e = ebegin; e != eend; ++e) {
				for (auto nr = e->begin(); nr != e->end(); ++nr) {
					for (auto nc = e->begin(); nc != e->end(); ++nc) {
						pattern.push_back({*nr, *nc});
					}
				}

//				for (auto nr = e->begin(), cbegin = e->begin(); nr != e->end(); ++nr, ++cbegin) {
//					for (auto nc = cbegin; nc != e->end(); ++nc) {
//						pattern.push_back(*nr <= *nc ? IJ{*nr, *nc} : IJ{*nc, *nr});
//					}
//				}
			}
			std::vector<eslocal> permutation(pattern.size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
				return pattern[i].row == pattern[j].row ? pattern[i].column < pattern[j].column : pattern[i].row < pattern[j].row;
			});

			DOFs.resize(pattern.size());
			for (size_t i = 1, nonzeros = 0; i < permutation.size(); i++) {
//				std::cout << "[" << pattern[permutation[i]].row + 1 << ", " << pattern[permutation[i]].column + 1 << "] = ";
				DOFs[permutation[i]] = pattern[permutation[i]] != pattern[permutation[i - 1]] ? ++nonzeros : nonzeros;
//				std::cout << DOFs[permutation[i]] << "\n";
			}

			_DOFsPermutation[d].swap(DOFs);
		}
	}

//	std::cout << *_domainNodes << "\n";
	std::cout << _DOFsPermutation.front();
}

void PhysicsInVectors::computeValues()
{
	size_t threads = environment->OMP_NUM_THREADS;

	eslocal eindex = 0, nindex = 0, vindex = 0;
//	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		DenseMatrix Ke, fe;
		for (eslocal d = _mesh.elements->domainDistribution[t]; d != _mesh.elements->domainDistribution[t + 1]; ++d) {
			std::fill(_instance.K[d].CSR_V_values.begin(), _instance.K[d].CSR_V_values.end(), 0);
			for (eslocal e = _mesh.elements->elementsDistribution[d]; e < _mesh.elements->elementsDistribution[d + 1]; ++e) {
				eslocal nsize = processElement(eindex++, nindex, Ke, fe);

//				for (auto r = 0, cbegin = 0; r < nsize; ++r, ++cbegin) {
//					for (auto c = cbegin; c < nsize; ++c, ++vindex) {
//						_instance.K[d].CSR_V_values[_DOFsPermutation[d][vindex]] += Ke(r, c);
//					}
//				}
				for (auto r = 0; r < nsize; ++r) {
					for (auto c = 0; c < nsize; ++c, ++vindex) {
						_instance.K[d].CSR_V_values[_DOFsPermutation[d][vindex]] += Ke(r, c);
					}
				}
				nindex += nsize;
			}
		}
	}
}


void PhysicsInVectors::buildGlobalCSRPattern()
{
	std::vector<IJ> KPattern;
	std::vector<eslocal> RHSPattern;

	_globalIndices.resize(_mesh.nodes->size);
	for (auto it = _mesh.nodes->pintervals.begin(); it != _mesh.nodes->pintervals.end(); ++it) {
		std::iota(_globalIndices.begin() + it->begin, _globalIndices.begin() + it->end, it->globalOffset);
	}

	for (auto e = _mesh.elements->procNodes->begin(); e != _mesh.elements->procNodes->end(); ++e) {
		for (auto nr = e->begin(); nr != e->end(); ++nr) {
			RHSPattern.push_back(_globalIndices[*nr]);
			for (auto nc = e->begin(); nc != e->end(); ++nc) {
				KPattern.push_back({ _globalIndices[*nr], _globalIndices[*nc] });
			}
		}
	}
	std::vector<eslocal> pK(KPattern.size());
	std::iota(pK.begin(), pK.end(), 0);
	std::sort(pK.begin(), pK.end(), [&] (eslocal i, eslocal j) {
		return KPattern[i].row == KPattern[j].row ? KPattern[i].column < KPattern[j].column : KPattern[i].row < KPattern[j].row;
	});

	std::vector<eslocal> pRHS(RHSPattern.size());
	std::iota(pRHS.begin(), pRHS.end(), 0);
	std::sort(pRHS.begin(), pRHS.end(), [&] (eslocal i, eslocal j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	std::vector<std::vector<IJ> > sKBuffer(_mesh.neighbours.size()), rKBuffer(_mesh.neighbours.size());
	std::vector<std::vector<eslocal> > sRHSBuffer(_mesh.neighbours.size()), rRHSBuffer(_mesh.neighbours.size());

	auto iK = pK.begin();
	auto iRHS = pRHS.begin();
	size_t n = 0;
	for (auto it = _mesh.nodes->pintervals.begin(); it != _mesh.nodes->pintervals.end() && it->sourceProcess < environment->MPIrank; ++it) {
		if (_mesh.neighbours[n] < it->sourceProcess) {
			++n;
		}
		while (RHSPattern[*iRHS] < it->end) {
			sRHSBuffer[0].push_back(RHSPattern[*iRHS++]);
		}
		while (KPattern[*iK].row < it->end) {
			sKBuffer[n].push_back(KPattern[*iK++]);
		}
		_instance.K.front().haloRows = it->end;
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
		Esutils::sortAndRemoveDuplicity(rRHSBuffer[i]);
		_neighRHSSize.push_back(rRHSBuffer[i].size());
	}

	_localK = pK.size();
	_localRHS = pRHS.size();
	pK.resize(KPattern.size());
	pRHS.resize(RHSPattern.size());
	std::iota(pK.begin() + _localK, pK.end(), _localK);
	std::iota(pRHS.begin() + _localRHS, pRHS.end(), _localRHS);
	std::sort(pK.begin(), pK.end(), [&] (eslocal i, eslocal j) {
		return KPattern[i].row == KPattern[j].row ? KPattern[i].column < KPattern[j].column : KPattern[i].row < KPattern[j].row;
	});
	std::sort(pRHS.begin(), pRHS.end(), [&] (eslocal i, eslocal j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	_pK.resize(KPattern.size());
	_pRHS.resize(RHSPattern.size());
	_instance.K.front().rows = _mesh.nodes->size;
	_instance.K.front().cols = _mesh.nodes->uniqueTotalSize;
	_instance.K.front().CSR_I_row_indices.resize(1);
	_instance.K.front().CSR_J_col_indices.resize(1);
	_instance.K.front().CSR_I_row_indices.reserve(_mesh.nodes->size + 1);
	_instance.K.front().CSR_J_col_indices.reserve(KPattern.size());

	_instance.K.front().CSR_I_row_indices.front() = 1;
	_instance.K.front().CSR_J_col_indices.front() = KPattern[pK.front()].column + 1;
	_pK[pK.front()] = 0;
	_pRHS[pRHS.front()] = 0;
	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
			++nonzeros;
			_instance.K.front().CSR_J_col_indices.push_back(KPattern[pK[i]].column + 1);
			if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
				_instance.K.front().CSR_I_row_indices.push_back(nonzeros + 1);
			}
		}
		_pK[pK[i]] = nonzeros;
	}
	for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
			++nonzeros;
		}
		_pRHS[pRHS[i]] = nonzeros;
	}
	_instance.K.front().CSR_I_row_indices.push_back(_instance.K.front().CSR_J_col_indices.size() + 1);
	_instance.K.front().CSR_V_values.resize(_instance.K.front().CSR_J_col_indices.size());
	_instance.K.front().nnz =_instance.K.front().CSR_J_col_indices.size() + 1;

	_instance.f.front().resize(_instance.K.front().rows);
	_instance.primalSolution.front().resize(_instance.K.front().rows);


//	std::cout << _globalDOFsPermutation;
}

void PhysicsInVectors::computeGlobalValues()
{
	eslocal nindex = 0, KIndex = 0, RHSIndex = 0;
	DenseMatrix Ke, fe;
	std::fill(_instance.K.front().CSR_V_values.begin(), _instance.K.front().CSR_V_values.end(), 0);
	std::fill(_instance.f.front().begin(), _instance.f.front().end(), 0);
	for (eslocal e = 0; e < _mesh.elements->size; ++e) {
		eslocal nsize = processElement(e, nindex, Ke, fe);
		for (auto r = 0; r < nsize; ++r, ++RHSIndex) {
			_instance.f.front()[_pRHS[RHSIndex]] += fe(r, 0);
			for (auto c = 0; c < nsize; ++c, ++KIndex) {
				_instance.K.front().CSR_V_values[_pK[KIndex]] += Ke(r, c);
			}
		}
		nindex += nsize;
	}
}

void PhysicsInVectors::synchronize()
{
//	std::cout << _instance.K.front();
//	std::cout << _instance.f.front();

	auto &COL = _instance.K.front().CSR_J_col_indices;
	auto &VAL = _instance.K.front().CSR_V_values;
	auto &RHS = _instance.f.front();

	std::vector<eslocal> ROW;
	for (size_t r = 0; r < _instance.K.front().CSR_I_row_indices.size() - 1; r++) {
		ROW.insert(ROW.end(), _instance.K.front().CSR_I_row_indices[r + 1] - _instance.K.front().CSR_I_row_indices[r], _globalIndices[r] + 1);
	}

//	Communication::serialize([&] () {
//		printf(" -- %d -- \n", environment->MPIrank);
//		for (eslocal r = 0, i = 0; r < _mesh.nodes->uniqueTotalSize; r++) {
//			for (eslocal c = 0; c < _mesh.nodes->uniqueTotalSize; c++) {
//				if (ROW[i] == r + 1 && COL[i] == c + 1) {
//					if (VAL[i] > -0.00001) {
//						if (VAL[i] > 10) {
//							printf(" %3.1f ", VAL[i++]);
//						} else {
//							printf(" %3.2f ", VAL[i++]);
//						}
//					} else {
//						printf("%3.2f ", VAL[i++]);
//					}
//				} else {
//					printf("      ");
//				}
//			}
//			printf("\n");
//		}
//	});

	std::vector<std::vector<double> > sBuffer(_mesh.neighbours.size()), rBuffer(_mesh.neighbours.size());

	size_t n = 0;
	for (auto it = _mesh.nodes->pintervals.begin(); it != _mesh.nodes->pintervals.end() && it->sourceProcess < environment->MPIrank; ++it) {
		if (_mesh.neighbours[n] < it->sourceProcess) {
			++n;
		}
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

//	if (environment->MPIrank == 0) {
//		std::cout << rBuffer.front();
//	}

	size_t KIndex = _localK, RHSIndex = _localRHS;
	for (size_t i = 0, j = 0; i < rBuffer.size(); ++i, j = 0) {
		for (; j < _neighRHSSize[i]; ++j, ++RHSIndex) {
			_instance.f.front()[_pRHS[RHSIndex]] += rBuffer[i][j];
		}
		for (; j < rBuffer[i].size(); ++j, ++KIndex) {
			_instance.K.front().CSR_V_values[_pK[KIndex]] += rBuffer[i][j];
		}
	}

//	Communication::serialize([&] () {
//		printf(" XX %d XX \n", environment->MPIrank);
//		for (eslocal r = 0, i = 0; r < _mesh.nodes->uniqueTotalSize; r++) {
//			for (eslocal c = 0; c < _mesh.nodes->uniqueTotalSize; c++) {
//				if (ROW[i] == r + 1 && COL[i] == c + 1) {
//					if (VAL[i] > -0.00001) {
//						if (VAL[i] > 10) {
//							printf(" %3.1f ", VAL[i++]);
//						} else {
//							printf(" %3.2f ", VAL[i++]);
//						}
//					} else {
//						printf("%3.2f ", VAL[i++]);
//					}
//				} else {
//					printf("      ");
//				}
//			}
//			printf("\n");
//		}
//	});

//	std::cout << _instance.K.front();
//	std::cout << _instance.f.front();
}







