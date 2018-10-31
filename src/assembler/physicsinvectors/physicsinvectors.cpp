
#include "physicsinvectors.h"

#include "../instance.h"

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
: _name(name), _mesh(mesh), _instance(instance), _step(step), _configuration(configuration), _DOFs{0, 0, 0, 1, 1}, _domainNodes(NULL)
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

void PhysicsInVectors::buildGlobalCSRPattern()
{
	std::vector<IJ> pattern;

	for (auto e = _mesh.elements->procNodes->begin(); e != _mesh.elements->procNodes->end(); ++e) {
		for (auto nr = e->begin(); nr != e->end(); ++nr) {
			for (auto nc = e->begin(); nc != e->end(); ++nc) {
				pattern.push_back({*nr, *nc});
			}
		}
	}
	std::vector<eslocal> permutation(pattern.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		return pattern[i].row == pattern[j].row ? pattern[i].column < pattern[j].column : pattern[i].row < pattern[j].row;
	});

	_globalDOFsPermutation.resize(pattern.size());
	_instance.K.front().rows = _mesh.nodes->size;
	_instance.K.front().cols = _mesh.nodes->size;
	_instance.K.front().CSR_I_row_indices.resize(1);
	_instance.K.front().CSR_J_col_indices.resize(1);
	_instance.K.front().CSR_I_row_indices.reserve(_mesh.nodes->size + 1);
	_instance.K.front().CSR_J_col_indices.reserve(pattern.size());

	_instance.K.front().CSR_I_row_indices.front() = 1;
	_instance.K.front().CSR_J_col_indices.front() = pattern[permutation.front()].column + 1;
	for (size_t i = 1, nonzeros = 0; i < permutation.size(); i++) {
		if (pattern[permutation[i]] != pattern[permutation[i - 1]]) {
			++nonzeros;
			_instance.K.front().CSR_J_col_indices.push_back(pattern[permutation[i]].column + 1);
			if (pattern[permutation[i - 1]].row != pattern[permutation[i]].row) {
				_instance.K.front().CSR_I_row_indices.push_back(nonzeros + 1);
			}
		}
		_globalDOFsPermutation[permutation[i]] = nonzeros;
	}
	_instance.K.front().CSR_I_row_indices.push_back(_instance.K.front().CSR_J_col_indices.size() + 1);
	_instance.K.front().CSR_V_values.resize(_instance.K.front().CSR_J_col_indices.size());
	_instance.K.front().nnz =_instance.K.front().CSR_J_col_indices.size() + 1;

	_instance.f.front().resize(_instance.K.front().rows);
	_instance.primalSolution.front().resize(_instance.K.front().rows);


//	std::cout << _globalDOFsPermutation;
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

void PhysicsInVectors::computeGlobalValues()
{
	eslocal nindex = 0, vindex = 0;
	DenseMatrix Ke, fe;
	std::fill(_instance.K.front().CSR_V_values.begin(), _instance.K.front().CSR_V_values.end(), 0);
	std::fill(_instance.f.front().begin(), _instance.f.front().end(), 0);
	for (eslocal e = 0; e < _mesh.elements->size; ++e) {
		eslocal nsize = processElement(e, nindex, Ke, fe);
		for (auto r = 0; r < nsize; ++r) {
			for (auto c = 0; c < nsize; ++c, ++vindex) {
				_instance.K.front().CSR_V_values[_globalDOFsPermutation[vindex]] += Ke(r, c);
			}
		}
		nindex += nsize;
	}
}

