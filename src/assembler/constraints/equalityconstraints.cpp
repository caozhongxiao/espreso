
#include "equalityconstraints.h"

#include "../instance.h"
#include "../step.h"

#include "../../basis/evaluators/evaluator.h"
#include "../../basis/containers/serializededata.h"

#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../mesh/mesh.h"
#include "../../mesh/store/boundaryregionstore.h"

#include "../../mesh/store/domainstore.h"
#include "../../config/ecf/environment.h"

#include <numeric>
#include <algorithm>


using namespace espreso;

EqualityConstraints::EqualityConstraints(Instance &instance, Mesh &mesh, const std::vector<BoundaryRegionStore*> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling)
: _instance(instance), _mesh(mesh)
{
	for (eslocal d = 0; d < _mesh._domains->size; d++) {
		_instance.B1[d].cols = _instance.domainDOFCount[d];
		_instance.B0[d].cols = _instance.domainDOFCount[d];
	}

	update(dirichlet, DOFs, withRedundantMultiplier, withScaling);
}

void EqualityConstraints::update(const std::vector<BoundaryRegionStore*> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling)
{
	_dirichlet = dirichlet;
	_DOFs = DOFs;
	_withRedundantMultipliers = withRedundantMultiplier;
	_withScaling = withScaling;

	_mergedDirichletIndices.clear();
	for (size_t i = 0; i < _mesh._domains->nodesIntervals.size(); ++i) {
		std::vector<eslocal> uniqueNodes;
		std::vector<double> values;
		for (size_t r = 0; r < dirichlet.size(); r++) {
			auto begin = dirichlet[r]->nodes->datatarray().begin() + dirichlet[r]->nodesIntervals[i].begin;
			auto end = dirichlet[r]->nodes->datatarray().begin() + dirichlet[r]->nodesIntervals[i].end;
			uniqueNodes.insert(uniqueNodes.end(), begin, end);
			for (auto n = begin; n != end; ++n) {
				// TODO: MESH
				values.push_back(r * 375.15);
			}
		}
		std::vector<eslocal> permutation(uniqueNodes.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return uniqueNodes[i] < uniqueNodes[j]; });

		_mergedDirichletIndices.push_back({});
		_mergedDirichletValues.push_back({});
		for (size_t p = 0; p < permutation.size(); ++p) {
			_mergedDirichletIndices.back().push_back(uniqueNodes[permutation[p]]);
			_mergedDirichletValues.back().push_back(values[permutation[p]]);
		}

		size_t unique = 0;
		for (size_t n = 1; n < _mergedDirichletIndices.back().size(); ++n) {
			if (_mergedDirichletIndices.back()[unique] != _mergedDirichletIndices.back()[n]) {
				_mergedDirichletIndices.back()[++unique] = _mergedDirichletIndices.back()[n];
				_mergedDirichletValues.back()[unique] = _mergedDirichletValues.back()[n];
			} else {
				if (_mergedDirichletValues.back()[unique] != _mergedDirichletValues.back()[n]) {
					ESINFO(ERROR) << "Multiple dirichlet values for a node.";
				}
			}
		}
		_mergedDirichletIndices.back().resize(std::min(unique + 1, _mergedDirichletIndices.back().size()));
		_mergedDirichletValues.back().resize(std::min(unique + 1, _mergedDirichletValues.back().size()));
	}

	auto getDirichletSize = [&] (eslocal i) -> eslocal {
		if (_withRedundantMultipliers) {
			return (eslocal)(_mesh._domains->nodesIntervals[i].ndomains * _mergedDirichletIndices[i].size());
		} else {
			return (eslocal)_mergedDirichletIndices[i].size();
		}
	};

	auto setDirichletOffset = [&] (eslocal i, esglobal offset) {
		for (size_t r = 0; r < dirichlet.size(); r++) {
			dirichlet[r]->nodesIntervals[i].LMOffset = offset;
			for (eslocal d = 0; d < _mesh._domains->size; d++) {
				eslocal ii = std::min(i, (eslocal)dirichlet[r]->domainNodesIntervals[d].size() - 1);
				while (ii >= 0 && dirichlet[r]->nodesIntervals[i].clusterOffset != dirichlet[r]->domainNodesIntervals[d][ii].clusterOffset) { --ii; }
				if (ii >= 0) {
					dirichlet[r]->domainNodesIntervals[d][ii].LMOffset = offset + _mergedDirichletIndices[i].size() * dirichlet[r]->domainNodesIntervals[d][ii].globalDomainOffset;
				}
			}
		}
	};

	_dirichletSize = _mesh.computeIntervalsOffsets(_mesh._domains->nodesIntervals, getDirichletSize, setDirichletOffset);

	auto getGluingSize = [&] (eslocal i) -> eslocal {
		if (_mesh._domains->nodesIntervals[i].ndomains == 1) {
			return 0;
		}
		eslocal ndomains = _mesh._domains->nodesIntervals[i].ndomains;
		eslocal size = _mesh._domains->nodesIntervals[i].end - _mesh._domains->nodesIntervals[i].begin;
		size -= _mergedDirichletIndices[i].size();
		if (_withRedundantMultipliers) {
			return size * (ndomains * (ndomains - 1) / 2);
		} else {
			return size * (ndomains - 1);
		}

	};

	auto setGluingOffset = [&] (eslocal i, esglobal offset) {
		if (_mesh._domains->nodesIntervals[i].ndomains > 1) {
			_mesh._domains->nodesIntervals[i].LMOffset = _dirichletSize + offset;
			for (size_t d = 0; d < _mesh._domains->size; d++) {
				eslocal ii = std::min(i, (eslocal)_mesh._domains->domainNodesIntervals[d].size() - 1);
				while (ii >= 0 && _mesh._domains->nodesIntervals[i].clusterOffset != _mesh._domains->domainNodesIntervals[d][ii].clusterOffset) { --ii; }
				if (ii >= 0) {
					_mesh._domains->domainNodesIntervals[d][ii].LMOffset = _dirichletSize + offset;
				}
			}
		}
	};

	_gluingSize = _mesh.computeIntervalsOffsets(_mesh._domains->nodesIntervals, getGluingSize, setGluingOffset);
}

void EqualityConstraints::B1DirichletInsert(const Step &step)
{
	if (_dirichlet.size() == 0) {
		return;
	}
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh._domains->domainDistribution[t]; d < _mesh._domains->domainDistribution[t + 1]; d++) {
			for (size_t i = 0, ii = 0; i < _mesh._domains->domainNodesIntervals[d].size(); ++i, ++ii) {
				if (_withRedundantMultipliers || _mesh._domains->domainNodesIntervals[d][i].globalDomainOffset == 0) {

					while (_mesh._domains->nodesIntervals[ii].clusterOffset != _mesh._domains->domainNodesIntervals[d][i].clusterOffset) { ++ii; }

					_instance.B1[d].I_row_indices.insert(_instance.B1[d].I_row_indices.end(), _mergedDirichletIndices[ii].size(), 0);
					std::iota(_instance.B1[d].I_row_indices.end() - _mergedDirichletIndices[ii].size(), _instance.B1[d].I_row_indices.end(), _dirichlet.front()->domainNodesIntervals[d][i].LMOffset + 1);

					for (size_t n = 0; n < _mergedDirichletIndices[ii].size(); ++n) {
						_instance.B1[d].J_col_indices.push_back(_mesh._domains->domainNodesIntervals[d][i].domainOffset + _mergedDirichletIndices[ii][n] - _mesh._domains->domainNodesIntervals[d][i].clusterOffset + 1);
					}

					_instance.B1c[d].insert(_instance.B1c[d].end(), _mergedDirichletValues[ii].begin(), _mergedDirichletValues[ii].end());

					_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), _mergedDirichletIndices[ii].size(), 1);
					_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), _mergedDirichletIndices[ii].size(), 1);

					_instance.B1subdomainsMap[d].insert(_instance.B1subdomainsMap[d].end(), _mergedDirichletIndices[ii].size(), 0);
					std::iota(_instance.B1subdomainsMap[d].end() - _mergedDirichletIndices[ii].size(), _instance.B1subdomainsMap[d].end(), _dirichlet.front()->domainNodesIntervals[d][i].LMOffset);

					_instance.B1[d].nnz = _instance.B1[d].I_row_indices.size();
					_instance.B1[d].rows = _dirichletSize;

					_instance.LB[d].resize(_instance.B1[d].nnz, -std::numeric_limits<double>::infinity());
				}
			}
		}
	}

	for (size_t i = 0; i < _mesh._domains->nodesIntervals.size(); ++i) {
		for (size_t d = 0; d < _mesh._domains->size; d++) {

			eslocal ii = std::min(i, _mesh._domains->domainNodesIntervals[d].size() - 1);
			while (ii >= 0 && _mesh._domains->nodesIntervals[i].clusterOffset != _mesh._domains->domainNodesIntervals[d][ii].clusterOffset) { --ii; }

			if (ii >= 0 && (_withRedundantMultipliers || _mesh._domains->domainNodesIntervals[d][ii].globalDomainOffset == 0)) {
				for (size_t n = 0; n < _mergedDirichletIndices[i].size(); n++) {
					_instance.B1clustersMap.push_back({
						(esglobal)(_dirichlet.front()->domainNodesIntervals[d][ii].LMOffset + n),
						(esglobal)environment->MPIrank
					});
				}
			}
		}
	}
}

void EqualityConstraints::B1GlueElements(const Step &step)
{
	auto redundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter) {
		const EInterval &interval = _mesh._domains->domainNodesIntervals[d][i];
		esglobal LMindex = interval.LMOffset + LMcounter * (interval.ndomains * (interval.ndomains - 1) / 2);
		eslocal DOFindex = interval.domainOffset + (begin - interval.clusterOffset) + 1;
		for (eslocal n = begin; n != end; ++n, ++DOFindex, ++LMcounter) {
			for (eslocal lm = 0; lm < interval.globalDomainOffset; LMindex += interval.ndomains - ++lm) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.globalDomainOffset - lm);
			}
			for (eslocal lm = interval.globalDomainOffset + 1; lm < interval.ndomains; LMindex += interval.ndomains - lm++) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + lm - interval.globalDomainOffset);
			}
			_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), interval.globalDomainOffset, -1);
			_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), interval.ndomains - interval.globalDomainOffset - 1, 1);
			_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), interval.ndomains - 1, 1.0 / interval.ndomains);
		}
	};

	auto nonredundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter) {
		const EInterval &interval = _mesh._domains->domainNodesIntervals[d][i];
		esglobal LMindex = interval.LMOffset + LMcounter * (interval.ndomains - 1);
		eslocal DOFindex = interval.domainOffset + (begin - interval.clusterOffset) + 1;
		for (eslocal n = begin; n != end; ++n, ++DOFindex, ++LMcounter) {
			if (interval.globalDomainOffset) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.globalDomainOffset - 1);
				_instance.B1[d].V_values.push_back(-1);
				_instance.B1duplicity[d].push_back(1.0 / interval.ndomains);
			}
			if (interval.globalDomainOffset + 1 < interval.ndomains) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.globalDomainOffset);
				_instance.B1[d].V_values.push_back(1);
				_instance.B1duplicity[d].push_back(1.0 / interval.ndomains);
			}
		}
	};

	auto glue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter) {
		if (_withRedundantMultipliers) {
			redundantglue(d, i, begin, end, LMcounter);
		} else {
			nonredundantglue(d, i, begin, end, LMcounter);
		}
	};

	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh._domains->domainDistribution[t]; d < _mesh._domains->domainDistribution[t + 1]; d++) {
			for (size_t i = 0, ii = 0; i < _mesh._domains->domainNodesIntervals[d].size(); ++i, ++ii) {
				eslocal LMcounter = 0;
				eslocal current = _mesh._domains->domainNodesIntervals[d][i].begin;
				while (_mesh._domains->nodesIntervals[ii].clusterOffset != _mesh._domains->domainNodesIntervals[d][i].clusterOffset) { ++ii; }
				for (auto end = _mergedDirichletIndices[ii].begin(); end != _mergedDirichletIndices[ii].end(); current = *end + 1, ++end) {
					glue(d, i, current, *end, LMcounter);
				}
				glue(d, i, current, _mesh._domains->domainNodesIntervals[d][i].end, LMcounter);
			}
			_instance.B1[d].rows += _gluingSize;
			_instance.B1[d].nnz = _instance.B1[d].I_row_indices.size();
			_instance.B1c[d].resize(_instance.B1duplicity[d].size());
			for (size_t i = _instance.B1subdomainsMap[d].size(); i < _instance.B1[d].I_row_indices.size(); ++i) {
				_instance.B1subdomainsMap[d].push_back(_instance.B1[d].I_row_indices[i] - 1);
			}
		}
	}

	std::vector<esglobal> domainProcDistribution = _mesh._domains->gatherProcsDistribution();
	for (size_t i = 0; i < _mesh._domains->nodesIntervals.size(); ++i) {
		const EInterval &interval = _mesh._domains->nodesIntervals[i];
		esglobal LMindex = 0;
		std::vector<std::vector<esglobal> > nmap;

		for (auto n1 = interval.neighbors[0] == -1 ? interval.neighbors.begin() + 1 : interval.neighbors.begin(); n1 != interval.neighbors.end(); ++n1) {
			int fromRank = std::lower_bound(domainProcDistribution.begin(), domainProcDistribution.end(), *n1 + 1) - domainProcDistribution.begin() - 1;
			for (auto n2 = n1 + 1; n2 != interval.neighbors.end(); ++n2, ++LMindex) {
				int toRank = std::lower_bound(domainProcDistribution.begin(), domainProcDistribution.end(), *n2 + 1) - domainProcDistribution.begin() - 1;

				if (fromRank == environment->MPIrank) {
					nmap.push_back({ LMindex, environment->MPIrank });
					if (toRank != environment->MPIrank) {
						nmap.back().push_back(toRank);
					}
				} else if (toRank == environment->MPIrank) {
					nmap.push_back({ LMindex, environment->MPIrank });
					if (fromRank != environment->MPIrank) {
						nmap.back().push_back(fromRank);
					}
				}

				if (!_withRedundantMultipliers) {
					break;
				}
			}
			if (!_withRedundantMultipliers) {
				break;
			}
		}

		size_t isize = _mesh._domains->nodesIntervals[i].end - _mesh._domains->nodesIntervals[i].begin - _mergedDirichletIndices[i].size();
		LMindex = interval.LMOffset;
		eslocal noffset;
		if (_withRedundantMultipliers) {
			noffset = interval.ndomains * (interval.ndomains - 1) / 2;
		} else {
			noffset = interval.ndomains - 1;
		}
		for (size_t n = 0; n < isize; ++n, LMindex += noffset) {
			_instance.B1clustersMap.insert(_instance.B1clustersMap.end(), nmap.begin(), nmap.end());
			for (auto map = _instance.B1clustersMap.end() - nmap.size(); map != _instance.B1clustersMap.end(); ++map) {
				map->front() += LMindex;
			}
		}
	}
}

//void EqualityConstraints::insertElementGluingToB1(const Step &step, bool withRedundantMultiplier, bool withScaling)
//{
//	auto n2i = [ & ] (size_t neighbour) {
//		return std::lower_bound(_mesh._neighbours.begin(), _mesh._neighbours.end(), neighbour) - _mesh._neighbours.begin();
//	};
//
//	std::vector<esglobal> lambdasID = computeLambdasID(step, withRedundantMultiplier);
//
//	std::vector<eslocal> permutation(lambdasID.size());
//	std::iota(permutation.begin(), permutation.end(), 0);
//
//	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
//		return (lambdasID[i] < 0) ? false : (lambdasID[j] < 0) ? true : lambdasID[i] < lambdasID[j];
//	});
//
//	auto it = std::find_if(permutation.begin(), permutation.end(), [&] (eslocal i) { return lambdasID[i] == -1; });
//	permutation.resize(it - permutation.begin());
//
//
//	size_t threads = environment->OMP_NUM_THREADS;
//	std::vector<size_t> distribution = Esutils::getDistribution(threads, permutation.size());
//
//	// threads x domains x data
//	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(_instance.domains));
//	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(_instance.domains));
//	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(_instance.domains));
//	std::vector<std::vector<std::vector<double> > > dup(threads, std::vector<std::vector<double> >(_instance.domains));
//
//	std::vector<std::vector<std::vector<esglobal> > > cMap(threads);
//
//	auto findDomain = [&] (const OldElement *e, size_t d, size_t dof) -> eslocal {
//		auto &DOFIndices = e->DOFsIndices();
//		size_t c = 0, DOFsSize = DOFIndices.size() / e->domains().size();
//		for (size_t i = 0; i < e->domains().size(); i++) {
//			if (DOFIndices[i * DOFsSize + dof] != -1) {
//				if (d == c++) {
//					return e->domains()[i];
//				}
//			}
//		}
//		return 0;
//	};
//
//	std::vector<std::vector<double> > diagonals;
//	if (withScaling) {
//		diagonals.resize(permutation.size());
//		std::vector<std::vector<double> > D(_instance.domains);
//
//		#pragma omp parallel for
//		for  (size_t d = 0; d < _instance.domains; d++) {
//			D[d] = _instance.K[d].getDiagonal();
//		}
//
//		std::vector<std::vector<std::vector<double> > > sBuffer(threads, std::vector<std::vector<double> >(_mesh._neighbours.size()));
//		std::vector<std::vector<double> > rBuffer(_mesh._neighbours.size());
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
//
//				const OldElement *e = _gluedElements[permutation[i] / _gluedDOFs.size()];
//				size_t dof = _gluedDOFsMeshOffsets[permutation[i] % _gluedDOFs.size()];
//
//				for (auto c = e->clusters().begin(); c != e->clusters().end(); ++c) {
//					if (*c != environment->MPIrank) {
//						for (auto d = e->domains().begin(); d != e->domains().end(); d++) {
//							sBuffer[t][n2i(*c)].push_back(D[*d][e->DOFIndex(*d, dof)]);
//						}
//					}
//				}
//
//			}
//		}
//
//		for (size_t t = 1; t < threads; t++) {
//			for (size_t n = 0; n < sBuffer[t].size(); n++) {
//				sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//			}
//		}
//
//		if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh._neighbours)) {
//			ESINFO(ERROR) << "problem while exchange K diagonal in B1 scaling.";
//		}
//
//		std::vector<eslocal> nPointer(_mesh._neighbours.size());
//		for (size_t i = 0; i < diagonals.size(); i++) {
//
//			const OldElement *e = _gluedElements[permutation[i] / _gluedDOFs.size()];
//			size_t dof = _gluedDOFsMeshOffsets[permutation[i] % _gluedDOFs.size()];
//
//			for (auto c = e->clusters().begin(); c != e->clusters().end(); ++c) {
//				if (*c == environment->MPIrank) {
//					for (auto d = e->domains().begin(); d != e->domains().end(); d++) {
//						diagonals[i].push_back(D[*d][e->DOFIndex(*d, dof)]);
//					}
//				} else {
//					for (eslocal d = 0; d < e->DOFCounter(*c, dof); d++) {
//						diagonals[i].push_back(rBuffer[n2i(*c)][nPointer[n2i(*c)]++]);
//					}
//				}
//			}
//
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
//
//			const OldElement *e = _gluedElements[permutation[i] / _gluedDOFs.size()];
//			size_t dof = _gluedDOFsMeshOffsets[permutation[i] % _gluedDOFs.size()];
//			esglobal offset = 0;
//			double duplicity = 0;
//			if (withScaling) {
//				std::for_each(diagonals[i].begin(), diagonals[i].end(), [&] (double v) { duplicity += v; });
//			} else {
//				duplicity = e->numberOfGlobalDomainsWithDOF(dof);
//			}
//
//			eslocal diag1 = 0;
//			for (auto c1 = e->clusters().begin(); c1 != e->clusters().end(); ++c1) {
//				for (eslocal d1 = 0; d1 < e->DOFCounter(*c1, dof); d1++, diag1++) {
//
//					eslocal diag2 = diag1 + 1;
//					for (auto c2 = c1; c2 != e->clusters().end(); ++c2) {
//						for (eslocal d2 = (*c1 == *c2 ? d1 + 1 : 0); d2 < e->DOFCounter(*c2, dof); d2++, diag2++) {
//
//							if (*c1 == environment->MPIrank) {
//								eslocal d = findDomain(e, d1, dof);
//								rows[t][d].push_back(lambdasID[permutation[i]] + offset + 1);
//								cols[t][d].push_back(e->DOFIndex(d, dof) + 1);
//								vals[t][d].push_back(1);
//								if (withScaling) {
//									dup[t][d].push_back(diagonals[i][diag2] / duplicity);
//								} else {
//									dup[t][d].push_back(1 / duplicity);
//								}
//							}
//
//							if (*c2 == environment->MPIrank) {
//								eslocal d = findDomain(e, d2, dof);
//								rows[t][d].push_back(lambdasID[permutation[i]] + offset + 1);
//								cols[t][d].push_back(e->DOFIndex(d, dof) + 1);
//								vals[t][d].push_back(-1);
//								if (withScaling) {
//									dup[t][d].push_back(diagonals[i][diag1] / duplicity);
//								} else {
//									dup[t][d].push_back(1 / duplicity);
//								}
//							}
//
//							if (*c1 == environment->MPIrank || *c2 == environment->MPIrank) {
//								cMap[t].push_back({ lambdasID[permutation[i]] + offset });
//								if (*c1 == *c2) {
//									cMap[t].back().push_back(*c1);
//								} else if (*c1 == environment->MPIrank) {
//									cMap[t].back().push_back(*c1);
//									cMap[t].back().push_back(*c2);
//								} else {
//									cMap[t].back().push_back(*c2);
//									cMap[t].back().push_back(*c1);
//								}
//							}
//
//							offset++;
//						}
//					}
//					if (!withRedundantMultiplier) {
//						break;
//					}
//				}
//				if (!withRedundantMultiplier) {
//					break;
//				}
//			}
//
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t p = 0; p < _instance.domains; p++) {
//		for (size_t t = 0; t < threads; t++) {
//			_instance.B1[p].I_row_indices.insert(_instance.B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
//			_instance.B1[p].J_col_indices.insert(_instance.B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
//			_instance.B1[p].V_values.insert(_instance.B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
//			_instance.B1duplicity[p].insert(_instance.B1duplicity[p].end(), dup[t][p].begin(), dup[t][p].end());
//		}
//		_instance.B1[p].cols = _instance.domainDOFCount[p];
//		_instance.B1[p].nnz = _instance.B1[p].I_row_indices.size();
//		_instance.B1c[p].resize(_instance.B1[p].nnz, 0);
//		_instance.LB[p].resize(_instance.B1[p].nnz, -std::numeric_limits<double>::infinity());
//		for (eslocal r = _instance.B1subdomainsMap[p].size(); r < _instance.B1[p].nnz; r++) {
//			_instance.B1subdomainsMap[p].push_back(_instance.B1[p].I_row_indices[r] - 1);
//		}
//	}
//
//	for (size_t t = 0; t < threads; t++) {
//		_instance.B1clustersMap.insert(_instance.B1clustersMap.end(), cMap[t].begin(), cMap[t].end());
//	}
//
//	ESINFO(EXHAUSTIVE) << "Lambdas in B1: " << _instance.B1[0].rows;
//}

void EqualityConstraints::insertCornersGluingToB0()
{
//	const std::vector<OldElement*> corners;
//	// const std::vector<OldElement*> &corners = _mesh.corners();
//	if (!corners.size()) {
//		ESINFO(ERROR) << "Mesh not contains corners.";
//		return;
//	}
//
//	for (size_t d = 0; d < _instance.domains; d++) {
//		_instance.B0[d].cols = _instance.K[d].cols;
//	}
//
//	size_t lambdas = _instance.B0[0].rows;
//
//	for (size_t e = 0; e < corners.size(); e++) {
//		for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
//			if (corners[e]->numberOfLocalDomainsWithDOF(_gluedDOFsMeshOffsets[dof]) > 1) { // inner nodes are not glued
//
//				for (size_t d1 = 0, d2 = 1; d2 < corners[e]->domains().size(); d1++, d2++) {
//
//					_instance.B0[corners[e]->domains()[d1]].I_row_indices.push_back(lambdas + 1);
//					_instance.B0[corners[e]->domains()[d1]].J_col_indices.push_back(corners[e]->DOFIndex(corners[e]->domains()[d1], _gluedDOFsMeshOffsets[dof]) + 1);
//					_instance.B0[corners[e]->domains()[d1]].V_values.push_back(1);
//
//					_instance.B0[corners[e]->domains()[d2]].I_row_indices.push_back(lambdas + 1);
//					_instance.B0[corners[e]->domains()[d2]].J_col_indices.push_back(corners[e]->DOFIndex(corners[e]->domains()[d2], _gluedDOFsMeshOffsets[dof]) + 1);
//					_instance.B0[corners[e]->domains()[d2]].V_values.push_back(-1);
//
//					lambdas++;
//				}
//
//			}
//		}
//	}
//
//	#pragma omp parallel for
//	for  (size_t p = 0; p < _instance.domains; p++) {
//		_instance.B0[p].rows = lambdas;
//		_instance.B0[p].cols = _instance.domainDOFCount[p];
//		_instance.B0[p].nnz = _instance.B0[p].I_row_indices.size();
//
//		_instance.B0subdomainsMap[p].reserve(_instance.B0[p].nnz);
//		for (eslocal i = _instance.B0subdomainsMap[p].size(); i < _instance.B0[p].nnz; i++) {
//			_instance.B0subdomainsMap[p].push_back(_instance.B0[p].I_row_indices[i] - 1);
//		}
//	}

	// TODO: MESH
	// ESINFO(EXHAUSTIVE) << "Average number of lambdas in B0 is " << Info::averageValue(lambdas);
}

void EqualityConstraints::insertKernelsGluingToB0(const std::vector<SparseMatrix> &kernels)
{
//	std::vector<OldElement*> el(_gluedInterfaceElements);
//
//	std::sort(el.begin(), el.end(), [] (OldElement* e1, OldElement* e2) {
//		if (e1->domains().size() != e2->domains().size()) {
//			return e1->domains().size() < e2->domains().size();
//		}
//		return e1->domains() < e2->domains();
//	});
//
//	std::vector<size_t> part;
//	part.push_back(std::lower_bound(el.begin(), el.end(), 2, [] (OldElement *e, size_t size) { return e->domains().size() < size; }) - el.begin());
//	ESTEST(MANDATORY) << "There are not elements on the sub-domains interface." << ((_gluedInterfaceElements.size() - part[0]) ? TEST_PASSED : TEST_FAILED);
//	for (size_t i = part[0] + 1; i < el.size(); i++) {
//		if (i && el[i - 1]->domains() != el[i]->domains()) {
//			part.push_back(i);
//		}
//	}
//	part.push_back(el.size());
//
//	std::vector<eslocal> rowIndex;
//	std::vector<eslocal> clusterRowIndex(_instance.clustersMap.size(), 1);
//	for (size_t i = 0; i < part.size() - 1; i++) {
//		const std::vector<eslocal> &domains = el[part[i]]->domains();
//		if (_instance.clustersMap[domains[0]] == _instance.clustersMap[domains[1]]) {
//			eslocal master = kernels[domains[0]].cols > kernels[domains[1]].cols ? domains[0] : domains[1];
//			eslocal rows = kernels[master].cols > 0 ? kernels[master].cols : 1;
//			rowIndex.push_back(clusterRowIndex[_instance.clustersMap[domains[0]]]);
//			clusterRowIndex[_instance.clustersMap[domains[0]]] += rows;
//		} else {
//			rowIndex.push_back(-1);
//		}
//	}
//
//	#pragma omp parallel for
//	for  (size_t p = 0; p < _instance.domains; p++) {
//		for (size_t i = 0; i < part.size() - 1; i++) {
//			const std::vector<eslocal> &domains = el[part[i]]->domains();
//			if (_instance.clustersMap[domains[0]] != _instance.clustersMap[domains[1]]) {
//				continue;
//			}
//
//			int sign = domains[0] == (eslocal)p ? 1 : domains[1] == (eslocal)p ? -1 : 0;
//			if (sign == 0) {
//				continue;
//			}
//
//			std::vector<OldElement*> DOFsOnInterface;
//			for (size_t e = part[i]; e < part[i + 1]; e++) {
//				if (_interfaceElementContainsGluedDOFs) {
//					for (size_t n = 0; n < el[e]->DOFsIndices().size(); n++) {
//						DOFsOnInterface.push_back(_gluedElements[el[e]->DOFsIndices()[n]]);
//					}
//				} else {
//					for (size_t n = 0; n < el[e]->nodes(); n++) {
//						DOFsOnInterface.push_back(_gluedElements[el[e]->node(n)]);
//					}
//				}
//			}
//			std::sort(DOFsOnInterface.begin(), DOFsOnInterface.end());
//			Esutils::removeDuplicity(DOFsOnInterface);
//
//			eslocal master = kernels[domains[0]].cols > kernels[domains[1]].cols ? domains[0] : domains[1];
//			if (kernels[master].cols == 0) {
//				for (size_t n = 0; n < DOFsOnInterface.size(); n++) {
//					for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
//						_instance.B0[p].I_row_indices.push_back(rowIndex[i]);
//						_instance.B0[p].J_col_indices.push_back(DOFsOnInterface[n]->DOFIndex(p, _gluedDOFsMeshOffsets[dof]) + 1);
//						_instance.B0[p].V_values.push_back(sign);
//					}
//				}
//			} else {
//				for (eslocal col = 0; col < kernels[master].cols; col++) {
//					for (size_t n = 0; n < DOFsOnInterface.size(); n++) {
//						for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
//							_instance.B0[p].I_row_indices.push_back(rowIndex[i] + col);
//							_instance.B0[p].J_col_indices.push_back(DOFsOnInterface[n]->DOFIndex(p, _gluedDOFsMeshOffsets[dof]) + 1);
//							_instance.B0[p].V_values.push_back(sign * kernels[master].dense_values[kernels[master].rows * col + DOFsOnInterface[n]->DOFIndex(master, _gluedDOFsMeshOffsets[dof])]);
//						}
//					}
//				}
//			}
//		}
//
//
//		_instance.B0[p].rows = clusterRowIndex[_instance.clustersMap[p]] - 1;
//		_instance.B0[p].cols = _instance.domainDOFCount[p];
//		_instance.B0[p].nnz = _instance.B0[p].I_row_indices.size();
//		_instance.B0subdomainsMap[p].reserve(_instance.B0[p].nnz);
//		for (eslocal i = _instance.B0subdomainsMap[p].size(); i < _instance.B0[p].nnz; i++) {
//			_instance.B0subdomainsMap[p].push_back(_instance.B0[p].I_row_indices[i] - 1);
//		}
//	}
}

#ifndef HAVE_MORTAR

void EqualityConstraints::insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave)
{
	ESINFO(GLOBAL_ERROR) << "Link 'mortarc' library.";
}

#else

#include "mortar.h"

void EqualityConstraints::insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave)
{

	std::vector<int> rows;
	std::vector<int> columns;
	std::vector<double> values;

	std::vector<std::vector<int> > masterElements;
	std::vector<Point_3D> masterCoordinates;
	std::vector<std::vector<int> > slaveElements;
	std::vector<Point_3D> slaveCoordinates;
	std::vector<int> nodes;

//	Region *mregion = _mesh.region(master);
//	Region *sregion = _mesh.region(slave);
//
//	std::vector<eslocal> masterElementsSBuffer, masterElementsRBuffer, slaveElementsSBuffer, slaveElementsRBuffer;
//	std::vector<esglobal> masterCoordinateIndicesSBuffer, slaveCoordinateIndicesSBuffer, masterCoordinateIndicesRBuffer, slaveCoordinateIndicesRBuffer;
//	std::vector<Point_3D> masterCoordinateSBuffer, slaveCoordinateSBuffer, masterCoordinateRBuffer, slaveCoordinateRBuffer;
//
//	for (size_t e = 0; e < mregion->elements().size(); e++) {
//		masterElementsSBuffer.push_back(mregion->elements()[e]->nodes());
//		for (size_t n = 0; n < mregion->elements()[e]->nodes(); n++) {
//			masterElementsSBuffer.push_back(_mesh.coordinates().globalIndex(mregion->elements()[e]->node(n)));
//			nodes.push_back(mregion->elements()[e]->node(n));
//		}
//	}
//
//	std::sort(nodes.begin(), nodes.end());
//	Esutils::removeDuplicity(nodes);
//
//	for (size_t n = 0; n < nodes.size(); n++) {
//		masterCoordinateSBuffer.push_back(Point_3D());
//		masterCoordinateSBuffer.back().x = _mesh.coordinates()[nodes[n]].x;
//		masterCoordinateSBuffer.back().y = _mesh.coordinates()[nodes[n]].y;
//		masterCoordinateSBuffer.back().z = _mesh.coordinates()[nodes[n]].z;
//		masterCoordinateIndicesSBuffer.push_back(_mesh.coordinates().globalIndex(nodes[n]));
//	}
//
//	nodes.clear();
//	for (size_t e = 0; e < sregion->elements().size(); e++) {
//		slaveElementsSBuffer.push_back(sregion->elements()[e]->nodes());
//		for (size_t n = 0; n < sregion->elements()[e]->nodes(); n++) {
//			slaveElementsSBuffer.push_back(_mesh.coordinates().globalIndex(sregion->elements()[e]->node(n)));
//			nodes.push_back(sregion->elements()[e]->node(n));
//		}
//	}
//
//	std::sort(nodes.begin(), nodes.end());
//	Esutils::removeDuplicity(nodes);
//
//	for (size_t n = 0; n < nodes.size(); n++) {
//		slaveCoordinateSBuffer.push_back(Point_3D());
//		slaveCoordinateSBuffer.back().x = _mesh.coordinates()[nodes[n]].x;
//		slaveCoordinateSBuffer.back().y = _mesh.coordinates()[nodes[n]].y;
//		slaveCoordinateSBuffer.back().z = _mesh.coordinates()[nodes[n]].z;
//		slaveCoordinateIndicesSBuffer.push_back(_mesh.coordinates().globalIndex(nodes[n]));
//	}

//	Communication::gatherUnknownSize(masterElementsSBuffer, masterElementsRBuffer);
//	Communication::gatherUnknownSize(slaveElementsSBuffer, slaveElementsRBuffer);
//	Communication::gatherUnknownSize(masterCoordinateSBuffer, masterCoordinateRBuffer);
//	Communication::gatherUnknownSize(slaveCoordinateSBuffer, slaveCoordinateRBuffer);
//	Communication::gatherUnknownSize(masterCoordinateIndicesSBuffer, masterCoordinateIndicesRBuffer);
//	Communication::gatherUnknownSize(slaveCoordinateIndicesSBuffer, slaveCoordinateIndicesRBuffer);
//
//	for (size_t e = 0; e < masterElementsRBuffer.size(); e += masterElementsRBuffer[e] + 1) {
//		masterElements.push_back(std::vector<int>(masterElementsRBuffer.begin() + e + 1, masterElementsRBuffer.begin() + e + 1 + masterElementsRBuffer[e]));
//	}
//	for (size_t e = 0; e < slaveElementsRBuffer.size(); e += slaveElementsRBuffer[e] + 1) {
//		slaveElements.push_back(std::vector<int>(slaveElementsRBuffer.begin() + e + 1, slaveElementsRBuffer.begin() + e + 1 + slaveElementsRBuffer[e]));
//	}
//
//	std::vector<eslocal> masterPermutation(masterCoordinateIndicesRBuffer.size()), slavePermutation(slaveCoordinateIndicesRBuffer.size()), masterUnique, slaveUnique;
//	std::iota(masterPermutation.begin(), masterPermutation.end(), 0);
//	std::iota(slavePermutation.begin(), slavePermutation.end(), 0);
//	std::sort(masterPermutation.begin(), masterPermutation.end(), [&] (eslocal i, eslocal j) { return masterCoordinateIndicesRBuffer[i] < masterCoordinateIndicesRBuffer[j]; });
//	std::sort(slavePermutation.begin(), slavePermutation.end(), [&] (eslocal i, eslocal j) { return slaveCoordinateIndicesRBuffer[i] < slaveCoordinateIndicesRBuffer[j]; });
//
//	for (size_t i = 0; i < masterPermutation.size(); i++) {
//		if (i == 0 || masterPermutation[i - 1] != masterPermutation[i]) {
//			masterCoordinates.push_back(masterCoordinateRBuffer[masterPermutation[i]]);
//			masterUnique.push_back(masterCoordinateIndicesRBuffer[masterPermutation[i]]);
//		}
//	}
//
//	for (size_t i = 0; i < slavePermutation.size(); i++) {
//		if (i == 0 || slavePermutation[i - 1] != slavePermutation[i]) {
//			slaveCoordinates.push_back(slaveCoordinateRBuffer[slavePermutation[i]]);
//			slaveUnique.push_back(slaveCoordinateIndicesRBuffer[slavePermutation[i]]);
//		}
//	}
//
//	for (size_t e = 0; e < masterElements.size(); e++) {
//		for (size_t n = 0; n < masterElements[e].size(); n++) {
//			masterElements[e][n] = std::lower_bound(masterUnique.begin(), masterUnique.end(), masterElements[e][n]) - masterUnique.begin();
//		}
//	}
//	for (size_t e = 0; e < slaveElements.size(); e++) {
//		for (size_t n = 0; n < slaveElements[e].size(); n++) {
//			slaveElements[e][n] = std::lower_bound(slaveUnique.begin(), slaveUnique.end(), slaveElements[e][n]) - slaveUnique.begin();
//		}
//	}
//
//	std::vector<IJV> d, m, supports, normals;
//	if (mregion->elements().size()) {
//		computeMortarEqualityConstraints(d, m, supports, normals, masterElements, masterCoordinates, slaveElements, slaveCoordinates);
//	}
}

#endif

