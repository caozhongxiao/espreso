
#include "equalityconstraints.h"

#include "../instance.h"
#include "../step.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../mesh/mesh.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/boundaryregionstore.h"

#include "../../config/ecf/environment.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

EqualityConstraints::EqualityConstraints(Instance &instance, Mesh &mesh, const std::map<std::string, ECFExpression> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling)
: _instance(instance), _mesh(mesh)
{
	for (eslocal d = 0; d < _mesh.elements->ndomains; d++) {
		_instance.B1[d].cols = _instance.domainDOFCount[d];
		_instance.B0[d].cols = _instance.domainDOFCount[d];
	}

	update(dirichlet, DOFs, withRedundantMultiplier, withScaling);
}

eslocal EqualityConstraints::computeIntervalsOffsets(std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, eslocal)> setsize)
{
	std::vector<eslocal> domainsProcDistribution = _mesh.elements->gatherDomainsProcDistribution();

	auto n2i = [&] (size_t n) {
		return std::lower_bound(_mesh.neighbours.begin(), _mesh.neighbours.end(), n) - _mesh.neighbours.begin();
	};
	auto d2p = [&] (eslocal d) {
		return std::lower_bound(domainsProcDistribution.begin(), domainsProcDistribution.end(), d + 1) - domainsProcDistribution.begin() - 1;
	};

	std::vector<std::vector<esglobal> > sGlobalOffset(_mesh.neighbours.size()), rGlobalOffset(_mesh.neighbours.size());
	esglobal offset = 0;

	auto domains = _mesh.nodes->idomains->cbegin();
	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++domains) {
		if (_mesh.elements->firstDomain <= domains->front() && domains->front() < _mesh.elements->firstDomain + _mesh.elements->ndomains) {
			auto begin = std::lower_bound(domains->begin(), domains->end(), _mesh.elements->firstDomain + _mesh.elements->ndomains) - 1;
			for (auto n = begin; n != domains->end(); ++n) {
				auto next = std::lower_bound(begin, domains->end(), *n);
				if (d2p(*begin) < d2p(*next)) {
					sGlobalOffset[n2i(d2p(*next))].push_back(offset);
				}
				begin = next;
			}
			offset += getsize(i);
		} else {
			rGlobalOffset[n2i(d2p(domains->front()))].push_back(0);
		}
	}

	Communication::exscan(offset);
	for (size_t n = 0; n < sGlobalOffset.size(); n++) {
		for (size_t i = 0; i < sGlobalOffset[n].size(); i++) {
			sGlobalOffset[n][i] += offset;
		}
	}

	if (!Communication::receiveLowerKnownSize(sGlobalOffset, rGlobalOffset, _mesh.neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: compute nodes global offset.";
	}

	domains = _mesh.nodes->idomains->cbegin();
	std::vector<int> rDataOffset(_mesh.neighbours.size());
	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++domains) {
		if (_mesh.elements->firstDomain <= domains->front() && domains->front() < _mesh.elements->firstDomain + _mesh.elements->ndomains) {
			setsize(i, offset);
			offset += getsize(i);
		} else {
			setsize(i, rGlobalOffset[n2i(d2p(domains->front()))][rDataOffset[n2i(d2p(domains->front()))]++]);
		}
	}

	MPI_Bcast(&offset, sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, environment->MPICommunicator);

	return offset;
}

void EqualityConstraints::update(const std::map<std::string, ECFExpression> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling)
{
	_DOFs = DOFs;
	_withRedundantMultipliers = withRedundantMultiplier;
	_withScaling = withScaling;

	_mergedDirichletOffset.clear();
	_mergedDirichletIndices.clear();
	_mergedDirichletValues.clear();
	std::vector<BoundaryRegionStore*> dregions;
	std::vector<Evaluator*> devaluator;
	for (auto it = dirichlet.begin(); it != dirichlet.end(); ++it) {
		dregions.push_back(_mesh.bregion(it->first));
		devaluator.push_back(it->second.evaluator);
	}
	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i) {
		std::vector<eslocal> uniqueNodes;
		std::vector<double> values;
		for (size_t r = 0; r < dregions.size(); r++) {
			auto begin = dregions[r]->nodes->datatarray().begin() + dregions[r]->nintervals[i].begin;
			auto end = dregions[r]->nodes->datatarray().begin() + dregions[r]->nintervals[i].end;
			uniqueNodes.insert(uniqueNodes.end(), begin, end);
			values.insert(values.end(), begin, end);
			devaluator[r]->evaluate(end - begin, NULL, NULL, 0, values.data() + values.size() - (end - begin));
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
	_mergedDirichletOffset.resize(_mergedDirichletIndices.size());

	auto getDirichletSize = [&] (eslocal i) -> eslocal {
		if (_withRedundantMultipliers) {
			return (_mesh.nodes->idomains->cbegin() + i)->size() * _mergedDirichletIndices[i].size();
		} else {
			return _mergedDirichletIndices[i].size();
		}
	};

	auto setDirichletOffset = [&] (eslocal i, esglobal offset) {
		_mergedDirichletOffset[i] = offset;
	};

	_dirichletSize = this->computeIntervalsOffsets(getDirichletSize, setDirichletOffset);

	auto getGluingSize = [&] (eslocal i) -> eslocal {
		size_t ndomains = (_mesh.nodes->idomains->cbegin() + i)->size();
		if (ndomains == 1) {
			return 0;
		}
		eslocal size = _mesh.nodes->pintervals[i].end - _mesh.nodes->pintervals[i].begin;
		size -= _mergedDirichletIndices[i].size();
		if (_withRedundantMultipliers) {
			return size * (ndomains * (ndomains - 1) / 2);
		} else {
			return size * (ndomains - 1);
		}
	};

	auto setGluingOffset = [&] (eslocal i, esglobal offset) {
		_intervalGluingOffset[i] = _dirichletSize + offset;
	};

	_intervalGluingOffset.clear();
	_intervalGluingOffset.resize(_mesh.nodes->pintervals.size());
	_gluingSize = this->computeIntervalsOffsets(getGluingSize, setGluingOffset);

	#pragma omp parallel for
	for (eslocal d = 0; d < _mesh.elements->ndomains; ++d) {
		for (size_t i = 0; i < _mesh.nodes->gintervals[d].size(); i++) {
				_mesh.nodes->gintervals[d][i].LMOffset = _intervalGluingOffset[_mesh.nodes->gintervals[d][i].pindex];
		}
	}
}

void EqualityConstraints::B1DirichletInsert(const Step &step)
{
	if (_dirichletSize == 0) {
		return;
	}
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = _mesh.elements->domainDistribution[t]; d < _mesh.elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < _mesh.nodes->gintervals[d].size(); ++i) {
				const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
				if (_withRedundantMultipliers || interval.dindex == 0) {

					_instance.B1[d].I_row_indices.insert(_instance.B1[d].I_row_indices.end(), _mergedDirichletIndices[interval.pindex].size(), 0);
					std::iota(
							_instance.B1[d].I_row_indices.end() - _mergedDirichletIndices[interval.pindex].size(),
							_instance.B1[d].I_row_indices.end(),
							_mergedDirichletOffset[interval.pindex] + interval.dindex * _mergedDirichletIndices[interval.pindex].size() + 1);

					for (size_t n = 0; n < _mergedDirichletIndices[interval.pindex].size(); ++n) {
						_instance.B1[d].J_col_indices.push_back(interval.DOFOffset + _mergedDirichletIndices[interval.pindex][n] - interval.begin + 1);
					}

					_instance.B1c[d].insert(_instance.B1c[d].end(), _mergedDirichletValues[interval.pindex].begin(), _mergedDirichletValues[interval.pindex].end());

					_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), _mergedDirichletIndices[interval.pindex].size(), 1);
					_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), _mergedDirichletIndices[interval.pindex].size(), 1);

					_instance.B1subdomainsMap[d].insert(_instance.B1subdomainsMap[d].end(), _mergedDirichletIndices[interval.pindex].size(), 0);
					std::iota(
							_instance.B1subdomainsMap[d].end() - _mergedDirichletIndices[interval.pindex].size(),
							_instance.B1subdomainsMap[d].end(),
							_mergedDirichletOffset[interval.pindex] + interval.dindex * _mergedDirichletIndices[interval.pindex].size());

					_instance.B1[d].nnz = _instance.B1[d].I_row_indices.size();
					_instance.B1[d].rows = _dirichletSize;

					_instance.LB[d].resize(_instance.B1[d].nnz, -std::numeric_limits<double>::infinity());
				}
			}
		}
	}

	auto iranks = _mesh.nodes->iranks->cbegin();
	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++iranks) {
		eslocal d = 0;
		for (auto r = iranks->begin(); r != iranks->end(); ++r, ++d) {
			if (*r == environment->MPIrank) {
				for (size_t n = 0; n < _mergedDirichletIndices[i].size(); ++n) {
					_instance.B1clustersMap.push_back({
						(esglobal)(_mergedDirichletOffset[i] + d * _mergedDirichletIndices[i].size() + n),
						(esglobal)environment->MPIrank
					});
				}
			}
			if (!_withRedundantMultipliers) {
				break;
			}
		}
	}
}

void EqualityConstraints::B1GlueElements(const Step &step)
{
	auto redundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter) {
		const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
		esglobal LMindex = interval.LMOffset + LMcounter * (interval.ndomains * (interval.ndomains - 1) / 2);
		eslocal DOFindex = interval.DOFOffset + (begin - interval.begin) + 1;
		for (eslocal n = begin; n != end; ++n, ++DOFindex, ++LMcounter) {
			for (eslocal lm = 0; lm < interval.dindex; LMindex += interval.ndomains - ++lm) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.dindex - lm);
			}
			for (eslocal lm = interval.dindex + 1; lm < interval.ndomains; ++LMindex, ++lm) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + 1);
			}
			LMindex += ((interval.ndomains - interval.dindex - 1) * (interval.ndomains - interval.dindex - 2) / 2);
			_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), interval.dindex, -1);
			_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), interval.ndomains - interval.dindex - 1, 1);
			_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), interval.ndomains - 1, 1.0 / interval.ndomains);
		}
	};

	auto nonredundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter) {
		const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
		esglobal LMindex = interval.LMOffset + LMcounter * (interval.ndomains - 1);
		eslocal DOFindex = interval.DOFOffset + (begin - interval.begin) + 1;
		for (eslocal n = begin; n != end; ++n, ++DOFindex, ++LMcounter) {
			if (interval.dindex) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.dindex - 1);
				_instance.B1[d].V_values.push_back(-1);
				_instance.B1duplicity[d].push_back(1.0 / interval.ndomains);
			}
			if (interval.dindex + 1 < interval.ndomains) {
				_instance.B1[d].J_col_indices.push_back(DOFindex);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.dindex);
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
		for (size_t d = _mesh.elements->domainDistribution[t]; d < _mesh.elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < _mesh.nodes->gintervals[d].size(); ++i) {
				const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
				eslocal LMcounter = 0;
				eslocal current = interval.begin;
				for (auto end = _mergedDirichletIndices[interval.pindex].begin(); end != _mergedDirichletIndices[interval.pindex].end(); current = *end + 1, ++end) {
					glue(d, i, current, *end, LMcounter);
				}
				glue(d, i, current, interval.end, LMcounter);
			}
			_instance.B1[d].rows += _gluingSize;
			_instance.B1[d].nnz = _instance.B1[d].I_row_indices.size();
			_instance.B1c[d].resize(_instance.B1duplicity[d].size());
			for (size_t i = _instance.B1subdomainsMap[d].size(); i < _instance.B1[d].I_row_indices.size(); ++i) {
				_instance.B1subdomainsMap[d].push_back(_instance.B1[d].I_row_indices[i] - 1);
			}
			_instance.LB[d].resize(_instance.B1[d].nnz, -std::numeric_limits<double>::infinity());
		}
	}


	auto iranks = _mesh.nodes->iranks->cbegin();
	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++iranks) {
		esglobal LMindex = 0;
		std::vector<std::vector<esglobal> > nmap;

		for (auto r1 = iranks->begin(); r1 != iranks->end(); ++r1) {
			for (auto r2 = r1 + 1; r2 != iranks->end(); ++r2, ++LMindex) {

				if (*r1 == environment->MPIrank) {
					nmap.push_back({ LMindex, environment->MPIrank });
					if (*r2 != environment->MPIrank) {
						nmap.back().push_back(*r2);
					}
				} else if (*r2 == environment->MPIrank) {
					nmap.push_back({ LMindex, environment->MPIrank });
					if (*r1 != environment->MPIrank) {
						nmap.back().push_back(*r1);
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

		eslocal noffset;
		if (_withRedundantMultipliers) {
			noffset = iranks->size() * (iranks->size() - 1) / 2;
		} else {
			noffset = iranks->size() - 1;
		}

		eslocal isize = _mesh.nodes->pintervals[i].end - _mesh.nodes->pintervals[i].begin - _mergedDirichletIndices[i].size();
		for (eslocal n = 0, LMoffset = _intervalGluingOffset[i]; n < isize; ++n, LMoffset += noffset) {
			_instance.B1clustersMap.insert(_instance.B1clustersMap.end(), nmap.begin(), nmap.end());
			for (auto map = _instance.B1clustersMap.end() - nmap.size(); map != _instance.B1clustersMap.end(); ++map) {
				map->front() += LMoffset;
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

