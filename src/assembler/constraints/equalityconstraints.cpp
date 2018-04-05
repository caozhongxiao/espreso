
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
#include "../../mesh/store/fetidatastore.h"
#include "../../config/ecf/environment.h"
#include "../../config/ecf/physics/physics.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

EqualityConstraints::EqualityConstraints(Instance &instance, Mesh &mesh, const RegionMap<ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling)
: _instance(instance), _mesh(mesh), _DOFs(1)
{
	update(dirichlet, withRedundantMultiplier, withScaling);
}

EqualityConstraints::EqualityConstraints(Instance &instance, Mesh &mesh, const RegionMap<ECFExpressionOptionalVector> &dirichlet, int DOFs, bool withRedundantMultiplier, bool withScaling)
: _instance(instance), _mesh(mesh), _DOFs(DOFs)
{
	update(dirichlet, withRedundantMultiplier, withScaling);
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

	Communication::exscan(offset, MPITools::operations().sizeToOffsetsEslocal);
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

void EqualityConstraints::update(const RegionMap<ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling)
{
	_dirichlet.clear();
	_dirichlet.resize(1);
	_intervalDirichletNodes.clear();
	_intervalDirichletNodes.resize(1);
	_intervalDirichletNodes[0].resize(_mesh.nodes->pintervals.size());

	_domainDirichletSize.resize(_mesh.elements->ndomains);

	for (auto it = dirichlet.regions.begin(); it != dirichlet.regions.end(); ++it) {
		_dirichlet[0].push_back(std::make_pair(_mesh.bregion(it->first), it->second.evaluator));
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); i++) {
			_intervalDirichletNodes[0][i].insert(_intervalDirichletNodes[0][i].end(),
					_dirichlet[0].back().first->uniqueNodes->datatarray().begin() + _dirichlet[0].back().first->unintervals[i].begin,
					_dirichlet[0].back().first->uniqueNodes->datatarray().begin() + _dirichlet[0].back().first->unintervals[i].end);
		}
	}

	for (auto it = dirichlet.intersections.begin(); it != dirichlet.intersections.end(); ++it) {
		_dirichlet[0].push_back(std::make_pair(_mesh.ibregion(it->first), it->second.evaluator));
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); i++) {
			_intervalDirichletNodes[0][i].insert(_intervalDirichletNodes[0][i].end(),
					_dirichlet[0].back().first->uniqueNodes->datatarray().begin() + _dirichlet[0].back().first->unintervals[i].begin,
					_dirichlet[0].back().first->uniqueNodes->datatarray().begin() + _dirichlet[0].back().first->unintervals[i].end);
		}
	}

	update(withRedundantMultiplier, withScaling);
}

void EqualityConstraints::update(const RegionMap<ECFExpressionOptionalVector> &dirichlet, bool withRedundantMultiplier, bool withScaling)
{
	_dirichlet.clear();
	_dirichlet.resize(_DOFs);
	_intervalDirichletNodes.clear();
	_intervalDirichletNodes.resize(_DOFs);
	for (int d = 0; d < _DOFs; d++) {
		_intervalDirichletNodes[d].resize(_mesh.nodes->pintervals.size());
	}

	_domainDirichletSize.resize(_mesh.elements->ndomains);

	auto insert = [&] (int d) {
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); i++) {
			_intervalDirichletNodes[d][i].insert(_intervalDirichletNodes[d][i].end(),
					_dirichlet[d].back().first->uniqueNodes->datatarray().begin() + _dirichlet[d].back().first->unintervals[i].begin,
					_dirichlet[d].back().first->uniqueNodes->datatarray().begin() + _dirichlet[d].back().first->unintervals[i].end);
		}
	};

	for (auto it = dirichlet.regions.begin(); it != dirichlet.regions.end(); ++it) {
		if (it->second.all.value.size()) {
			for (int d = 0; d < _DOFs; d++) {
				_dirichlet[d].push_back(std::make_pair(_mesh.bregion(it->first), it->second.all.evaluator));
				insert(d);
			}
		} else {
			if (it->second.x.value.size()) {
				_dirichlet[0].push_back(std::make_pair(_mesh.bregion(it->first), it->second.x.evaluator));
				insert(0);
			}
			if (it->second.y.value.size()) {
				_dirichlet[1].push_back(std::make_pair(_mesh.bregion(it->first), it->second.y.evaluator));
				insert(1);
			}
			if (_DOFs == 3 && it->second.z.value.size()) {
				_dirichlet[2].push_back(std::make_pair(_mesh.bregion(it->first), it->second.z.evaluator));
				insert(2);
			}
		}
	}

	for (auto it = dirichlet.intersections.begin(); it != dirichlet.intersections.end(); ++it) {
		if (it->second.all.value.size()) {
			for (int d = 0; d < _DOFs; d++) {
				_dirichlet[d].push_back(std::make_pair(_mesh.ibregion(it->first), it->second.all.evaluator));
				insert(d);
			}
		} else {
			if (it->second.x.value.size()) {
				_dirichlet[0].push_back(std::make_pair(_mesh.ibregion(it->first), it->second.x.evaluator));
				insert(0);
			}
			if (it->second.y.value.size()) {
				_dirichlet[1].push_back(std::make_pair(_mesh.ibregion(it->first), it->second.y.evaluator));
				insert(1);
			}
			if (_DOFs == 3 && it->second.z.value.size()) {
				_dirichlet[2].push_back(std::make_pair(_mesh.ibregion(it->first), it->second.z.evaluator));
				insert(2);
			}
		}
	}

	update(withRedundantMultiplier, withScaling);
}

void EqualityConstraints::update(bool withRedundantMultiplier, bool withScaling)
{
	_withRedundantMultipliers = withRedundantMultiplier;
	_withScaling = withScaling;

	for (eslocal d = 0; d < _mesh.elements->ndomains; d++) {
		_instance.B1[d].cols = _instance.domainDOFCount[d];
		_instance.B0[d].cols = _instance.domainDOFCount[d];
	}

	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); i++) {
		for (size_t d = 0; d < _intervalDirichletNodes.size(); d++) {
			size_t presize = _intervalDirichletNodes[d][i].size();
			Esutils::sortAndRemoveDuplicity(_intervalDirichletNodes[d][i]);
			if (presize != _intervalDirichletNodes[d][i].size()) {
				ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: DIRICHLET intervals have to be disjunct.";
			}
		}
	}

	std::vector<eslocal> dsize(_mesh.nodes->pintervals.size());
	for (size_t dof = 0; dof < _dirichlet.size(); dof++) {
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); i++) {
			dsize[i] += _intervalDirichletNodes[dof][i].size();
		}
	}

	auto getDirichletSize = [&] (eslocal i) -> eslocal {
		if (_withRedundantMultipliers) {
			return (_mesh.nodes->idomains->cbegin() + i)->size() * dsize[i];
		} else {
			return dsize[i];
		}
	};

	auto setDirichletOffset = [&] (eslocal i, esglobal offset) {
		_intervalDirichletOffset[i] = offset;
	};

	_intervalDirichletOffset.clear();
	_intervalDirichletOffset.resize(_mesh.nodes->pintervals.size());
	_dirichletSize = this->computeIntervalsOffsets(getDirichletSize, setDirichletOffset);

	auto getGluingSize = [&] (eslocal i) -> eslocal {
		size_t ndomains = (_mesh.nodes->idomains->cbegin() + i)->size();
		if (ndomains == 1) {
			return 0;
		}
		eslocal size = _DOFs * (_mesh.nodes->pintervals[i].end - _mesh.nodes->pintervals[i].begin);
		if (_withRedundantMultipliers) {
			return (size - dsize[i]) * (ndomains * (ndomains - 1) / 2);
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
				eslocal LMoffset = 0;
				const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
				for (int dof = 0; dof < _DOFs; dof++) {
					if (_withRedundantMultipliers || interval.dindex == 0) {

						_instance.B1[d].I_row_indices.insert(_instance.B1[d].I_row_indices.end(), _intervalDirichletNodes[dof][interval.pindex].size(), 0);
						std::iota(
								_instance.B1[d].I_row_indices.end() - _intervalDirichletNodes[dof][interval.pindex].size(),
								_instance.B1[d].I_row_indices.end(),
								_intervalDirichletOffset[interval.pindex] + interval.dindex * _intervalDirichletNodes[dof][interval.pindex].size() + LMoffset + 1);

						for (size_t r = 0; r < _dirichlet[dof].size(); r++) {
							for (
									auto n = _dirichlet[dof][r].first->uniqueNodes->datatarray().cbegin() + _dirichlet[dof][r].first->unintervals[interval.pindex].begin;
									n != _dirichlet[dof][r].first->uniqueNodes->datatarray().cbegin() + _dirichlet[dof][r].first->unintervals[interval.pindex].end;
									++n) {
								_instance.B1[d].J_col_indices.push_back(_DOFs * (interval.DOFOffset + *n - interval.begin) + dof + 1);
							}
						}

						_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), _intervalDirichletNodes[dof][interval.pindex].size(), 1);
						_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), _intervalDirichletNodes[dof][interval.pindex].size(), 1);

						_instance.B1subdomainsMap[d].insert(_instance.B1subdomainsMap[d].end(), _intervalDirichletNodes[dof][interval.pindex].size(), 0);
						std::iota(
								_instance.B1subdomainsMap[d].end() - _intervalDirichletNodes[dof][interval.pindex].size(),
								_instance.B1subdomainsMap[d].end(),
								_intervalDirichletOffset[interval.pindex] + interval.dindex * _intervalDirichletNodes[dof][interval.pindex].size() + LMoffset);

						_instance.B1[d].nnz = _instance.B1[d].I_row_indices.size();
						_instance.B1[d].rows = _dirichletSize;

						_instance.LB[d].resize(_instance.B1[d].nnz, -std::numeric_limits<double>::infinity());

						_domainDirichletSize[d] = _instance.B1[d].nnz;

						LMoffset += interval.ndomains * _intervalDirichletNodes[dof][interval.pindex].size();
					}
				}
			}
			_instance.B1c[d].resize(_instance.B1[d].V_values.size());
		}
	}
	_instance.block[Instance::CONSTRAINT::DIRICHLET] = _dirichletSize;


	B1DirichletUpdate(step);

	auto iranks = _mesh.nodes->iranks->cbegin();
	for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++iranks) {
		eslocal d = 0;
		eslocal LMoffset = _intervalDirichletOffset[i];
		for (int dof = 0; dof < _DOFs; dof++) {
			for (auto r = iranks->begin(); r != iranks->end(); ++r, ++d) {
				if (*r == environment->MPIrank) {
					for (size_t n = 0; n < _intervalDirichletNodes[dof][i].size(); ++n) {
						_instance.B1clustersMap.push_back({
							(esglobal)(n + LMoffset),
							(esglobal)environment->MPIrank
						});
					}
				}
				LMoffset += _intervalDirichletNodes[dof][i].size();
				if (!_withRedundantMultipliers) {
					break;
				}
			}
		}
	}
}

void EqualityConstraints::B1DirichletUpdate(const Step &step)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = _mesh.elements->domainDistribution[t]; d < _mesh.elements->domainDistribution[t + 1]; d++) {
			eslocal offset = 0;
			for (size_t i = 0; i < _mesh.nodes->gintervals[d].size(); ++i) {
				const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
				for (int dof = 0; dof < _DOFs; dof++) {
					if (_withRedundantMultipliers || interval.dindex == 0) {
						for (size_t r = 0; r < _dirichlet[dof].size(); r++) {
							size_t isize = _dirichlet[dof][r].first->unintervals[interval.pindex].end - _dirichlet[dof][r].first->unintervals[interval.pindex].begin;
							if (_dirichlet[dof][r].second->isTemperatureDependent()) {
								ESINFO(ERROR) << "Dirichlet boundary condition cannot be dependent on TEMPERATURE.";
							}
							if (_dirichlet[dof][r].second->isCoordinateDependent()) {
								std::vector<Point> points;
								points.reserve(isize);
								for (eslocal n = _dirichlet[dof][r].first->unintervals[interval.pindex].begin; n < _dirichlet[dof][r].first->unintervals[interval.pindex].end; ++n) {
									points.push_back(_mesh.nodes->coordinates->datatarray()[_dirichlet[dof][r].first->uniqueNodes->datatarray()[n]]);
								}
								_dirichlet[dof][r].second->evaluate(isize, points.data(), NULL, step.currentTime, _instance.B1c[d].data() + offset);
							} else {
								_dirichlet[dof][r].second->evaluate(isize, NULL, NULL, step.currentTime, _instance.B1c[d].data() + offset);
							}
							if (step.internalForceReduction != 1) {
								for (eslocal i = 0; i < isize; i++) {
									_instance.B1c[d][offset + i] *= step.internalForceReduction;
								}
							}
							offset += isize;
						}
					}
				}
			}
		}
	}
}

void EqualityConstraints::B1GlueElements()
{
	auto redundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter, int dof) {
		const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
		esglobal LMindex = interval.LMOffset + LMcounter * (interval.ndomains * (interval.ndomains - 1) / 2);
		eslocal DOFindex = interval.DOFOffset + (begin - interval.begin);
		for (eslocal n = begin; n != end; ++n, ++DOFindex, ++LMcounter) {
			for (eslocal lm = 0; lm < interval.dindex; LMindex += interval.ndomains - ++lm) {
				_instance.B1[d].J_col_indices.push_back(_DOFs * DOFindex + dof + 1);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.dindex - lm);
			}
			for (eslocal lm = interval.dindex + 1; lm < interval.ndomains; ++LMindex, ++lm) {
				_instance.B1[d].J_col_indices.push_back(_DOFs * DOFindex + dof + 1);
				_instance.B1[d].I_row_indices.push_back(LMindex + 1);
			}
			LMindex += ((interval.ndomains - interval.dindex - 1) * (interval.ndomains - interval.dindex - 2) / 2);
			_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), interval.dindex, -1);
			_instance.B1[d].V_values.insert(_instance.B1[d].V_values.end(), interval.ndomains - interval.dindex - 1, 1);
		}
		_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), (end - begin) * (interval.ndomains - 1), 1.0 / interval.ndomains);
	};

	auto nonredundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter, int dof) {
		const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
		esglobal LMindex = interval.LMOffset + LMcounter * (interval.ndomains - 1);
		eslocal DOFindex = interval.DOFOffset + (begin - interval.begin);
		for (eslocal n = begin, lm = 0; n != end; ++n, ++lm, ++DOFindex, ++LMcounter) {
			if (interval.dindex) {
				_instance.B1[d].J_col_indices.push_back(_DOFs * DOFindex + dof + 1);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.dindex + lm);
				_instance.B1[d].V_values.push_back(-1);
			}
			if (interval.dindex + 1 < interval.ndomains) {
				_instance.B1[d].J_col_indices.push_back(_DOFs * DOFindex + dof + 1);
				_instance.B1[d].I_row_indices.push_back(LMindex + interval.dindex + 1 + lm);
				_instance.B1[d].V_values.push_back(1);
			}
		}
		if (interval.dindex) {
			_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), end - begin, 1.0 / interval.ndomains);
		}
		if (interval.dindex + 1 < interval.ndomains) {
			_instance.B1duplicity[d].insert(_instance.B1duplicity[d].end(), end - begin, 1.0 / interval.ndomains);
		}
	};

	auto glue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter, int dof) {
		if (_withRedundantMultipliers) {
			redundantglue(d, i, begin, end, LMcounter, dof);
		} else {
			nonredundantglue(d, i, begin, end, LMcounter, dof);
		}
	};

	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = _mesh.elements->domainDistribution[t]; d < _mesh.elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < _mesh.nodes->gintervals[d].size(); ++i) {
				eslocal LMcounter = 0;
				const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
				for (int dof = 0; dof < _DOFs; dof++) {
					if (_withRedundantMultipliers) {
						eslocal current = interval.begin;
						for (auto end = _intervalDirichletNodes[dof][interval.pindex].begin(); end != _intervalDirichletNodes[dof][interval.pindex].end(); current = *end + 1, ++end) {
							glue(d, i, current, *end, LMcounter, dof);
						}
						glue(d, i, current, interval.end, LMcounter, dof);
					} else {
						glue(d, i, _mesh.nodes->gintervals[d][i].begin, _mesh.nodes->gintervals[d][i].end, LMcounter, dof);
					}
				}
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
			for (auto r2 = r1 + 1; r2 != iranks->end(); ++r2) {
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
				++LMindex;
				if (!_withRedundantMultipliers) {
					break;
				}
			}
		}

		eslocal noffset;
		if (_withRedundantMultipliers) {
			noffset = iranks->size() * (iranks->size() - 1) / 2;
		} else {
			noffset = iranks->size() - 1;
		}

		eslocal LMoffset = _intervalGluingOffset[i];
		for (int dof = 0; dof < _DOFs; dof++) {
			eslocal isize = _mesh.nodes->pintervals[i].end - _mesh.nodes->pintervals[i].begin;
			if (_withRedundantMultipliers) {
				isize -= _intervalDirichletNodes[dof][i].size();
			}
			for (eslocal n = 0; n < isize; ++n, LMoffset += noffset) {
				_instance.B1clustersMap.insert(_instance.B1clustersMap.end(), nmap.begin(), nmap.end());
				for (auto map = _instance.B1clustersMap.end() - nmap.size(); map != _instance.B1clustersMap.end(); ++map) {
					map->front() += LMoffset;
				}
			}
		}
	}

	_instance.block[Instance::CONSTRAINT::EQUALITY_CONSTRAINTS] = _dirichletSize + _gluingSize;

	if (_withScaling) {
		B1DuplicityUpdate();
	}
}

void EqualityConstraints::B1DuplicityUpdate()
{
	if (!_withScaling || !_withRedundantMultipliers) {
		return;
	}

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours.begin(), _mesh.neighbours.end(), neighbour) - _mesh.neighbours.begin();
	};

	// interval x DOF0(d0, d1, ...), DOF1(d0, d1, ...)
	std::vector<std::vector<double> > diagonals(_mesh.nodes->pintervals.size());
	if (_withScaling) {
		std::vector<std::vector<double> > D(_mesh.elements->ndomains);
		#pragma omp parallel for
		for  (size_t d = 0; d < _instance.domains; d++) {
			D[d] = _instance.K[d].getDiagonal();
		}

		std::vector<std::vector<double> > sBuffer(_mesh.neighbours.size()), rBuffer(_mesh.neighbours.size());
		auto irank = _mesh.nodes->iranks->cbegin();
		auto idomain = _mesh.nodes->idomains->cbegin();
		std::vector<int> dindex(_mesh.elements->ndomains), ranks;
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++irank, ++idomain) {
			ranks.clear();
			ranks.insert(ranks.end(), irank->begin(), irank->end());
			Esutils::sortAndRemoveDuplicity(ranks);
			for (auto rank = ranks.begin(); rank != ranks.end(); ++rank) {
				if (*rank != environment->MPIrank) {
					eslocal target = n2i(*rank);
					for (auto domain = idomain->begin(); domain != idomain->end(); ++domain) {
						if (_mesh.elements->firstDomain <= *domain && *domain < _mesh.elements->firstDomain + _mesh.elements->ndomains) {
							eslocal d = *domain - _mesh.elements->firstDomain;
							auto begin = D[d].begin() + _DOFs * _mesh.nodes->dintervals[d][dindex[d]].DOFOffset;
							auto end = begin + _DOFs * (_mesh.nodes->dintervals[d][dindex[d]].end - _mesh.nodes->dintervals[d][dindex[d]].begin);
							sBuffer[target].insert(sBuffer[target].end(), begin, end);
						}
					}
				}
			}
			for (auto domain = idomain->begin(); domain != idomain->end(); ++domain) {
				if (_mesh.elements->firstDomain <= *domain && *domain < _mesh.elements->firstDomain + _mesh.elements->ndomains) {
					++dindex[*domain - _mesh.elements->firstDomain];
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _mesh.neighbours)) {
			ESINFO(ERROR) << "problem while exchange K diagonal in B1 scaling.";
		}

		irank = _mesh.nodes->iranks->cbegin();
		idomain = _mesh.nodes->idomains->cbegin();
		std::fill(dindex.begin(), dindex.end(), 0);
		std::vector<eslocal> rindex(_mesh.neighbours.size());
		for (size_t i = 0; i < _mesh.nodes->pintervals.size(); ++i, ++irank, ++idomain) {
			if (idomain->size() > 1) {
				eslocal isize = _mesh.nodes->pintervals[i].end - _mesh.nodes->pintervals[i].begin;
				diagonals[i].resize(_DOFs * isize * idomain->size());
				for (int r = 0; r < irank->size(); ++r) {
					if (irank->at(r) != environment->MPIrank) {
						eslocal target = n2i(irank->at(r));
						for (eslocal n = 0; n < isize; ++n) {
							for (int dof = 0; dof < _DOFs; dof++) {
								diagonals[i][_DOFs * n * idomain->size() + dof * idomain->size() + r] = rBuffer[target][rindex[target]++];
							}
						}
					} else {
						eslocal d = idomain->at(r) - _mesh.elements->firstDomain;
						for (eslocal n = 0; n < isize; ++n) {
							for (int dof = 0; dof < _DOFs; dof++) {
								diagonals[i][_DOFs * n * idomain->size() + dof * idomain->size() + r] = D[d][_DOFs * (_mesh.nodes->dintervals[d][dindex[d]].DOFOffset + n) + dof];
							}
						}
						++dindex[d];
					}
				}
			}
		}
	}

	std::vector<eslocal> dindex = _domainDirichletSize;

	auto redundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter, int dof) {
		const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
		double sum;
		for (eslocal n = begin; n != end; ++n) {
			sum = 0;
			for (eslocal i = 0; i < interval.ndomains; i++) {
				sum += diagonals[interval.pindex][_DOFs * (n - interval.begin) * interval.ndomains + dof * interval.ndomains + i];
			}
			for (eslocal i = 0; i < interval.dindex; ++i) {
				_instance.B1duplicity[d][dindex[d]++] = diagonals[interval.pindex][_DOFs * (n - interval.begin) * interval.ndomains + dof * interval.ndomains + i] / sum;
			}
			for (eslocal i = interval.dindex + 1; i < interval.ndomains; ++i) {
				_instance.B1duplicity[d][dindex[d]++] = diagonals[interval.pindex][_DOFs * (n - interval.begin) * interval.ndomains + dof * interval.ndomains + i] / sum;
			}
		}
	};

	auto nonredundantglue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter, int dof) {
		// CURRENT SCALING WORKS ONLY WITH REDUNDANT MULTIPLIERS
	};

	auto glue = [&] (eslocal d, size_t i, eslocal begin, eslocal end, eslocal &LMcounter, int dof) {
		if (_withRedundantMultipliers) {
			redundantglue(d, i, begin, end, LMcounter, dof);
		} else {
			nonredundantglue(d, i, begin, end, LMcounter, dof);
		}
	};

	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = _mesh.elements->domainDistribution[t]; d < _mesh.elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < _mesh.nodes->gintervals[d].size(); ++i) {
				if (_mesh.nodes->gintervals[d][i].ndomains > 1) {
					eslocal LMcounter = 0;
					for (int dof = 0; dof < _DOFs; dof++) {
						if (_withRedundantMultipliers) {
							const GluingInterval &interval = _mesh.nodes->gintervals[d][i];
							eslocal current = interval.begin;
							for (auto end = _intervalDirichletNodes[dof][interval.pindex].begin(); end != _intervalDirichletNodes[dof][interval.pindex].end(); current = *end + 1, ++end) {
								glue(d, i, current, *end, LMcounter, dof);
							}
							glue(d, i, current, interval.end, LMcounter, dof);
						} else {
							glue(d, i, _mesh.nodes->gintervals[d][i].begin, _mesh.nodes->gintervals[d][i].end, LMcounter, dof);
						}
					}
				}
			}
		}
	}
}

void EqualityConstraints::B0Kernels(const std::vector<SparseMatrix> &kernels)
{
	std::vector<eslocal> rowIndex(_mesh.FETIData->inodesDomains.size());
	std::vector<eslocal> rCounters(*std::max_element(_mesh.elements->clusters.begin(), _mesh.elements->clusters.end()) + 1);

	for (size_t i = 0; i < _mesh.FETIData->inodesDomains.size(); i++) {
		eslocal domain = _mesh.FETIData->inodesDomains[i].first;
		eslocal ndomain = _mesh.FETIData->inodesDomains[i].second;
		eslocal cluster = _mesh.elements->clusters[domain];
		rowIndex[i] = rCounters[cluster];
		rCounters[cluster] += std::max(kernels[domain].cols, std::max(kernels[ndomain].cols, (eslocal)1));
	}

	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &nodes = _mesh.FETIData->interfaceNodes->datatarray();
		int sign, cols, master;
		for (eslocal d = _mesh.elements->domainDistribution[t]; d < _mesh.elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < _mesh.FETIData->inodesDomains.size(); i++) {
				sign = 0;
				if (_mesh.FETIData->inodesDomains[i].first == d) {
					sign = 1;
				}
				if (_mesh.FETIData->inodesDomains[i].second == d) {
					sign = -1;
				}

				if (sign != 0) {
					master = _mesh.FETIData->inodesDomains[i].first;
					if (kernels[master].cols < kernels[_mesh.FETIData->inodesDomains[i].second].cols) {
						master = _mesh.FETIData->inodesDomains[i].second;
					}
					cols = kernels[master].cols;
					if (cols) {
						for (eslocal c = 0; c < cols; c++) {
							auto dit = _mesh.nodes->dintervals[d].begin();
							auto masterit = _mesh.nodes->dintervals[master].begin();
							for (size_t n = _mesh.FETIData->inodesDistribution[i]; n < _mesh.FETIData->inodesDistribution[i + 1]; n++) {
								while (dit->end < nodes[n]) {
									++dit;
								}
								while (masterit->end < nodes[n]) {
									++masterit;
								}
								for (int dof = 0; dof < _DOFs; dof++) {
									_instance.B0[d].I_row_indices.push_back(rowIndex[i] + c + 1);
									_instance.B0[d].J_col_indices.push_back(_DOFs * (dit->DOFOffset + nodes[n] - dit->begin) + dof + 1);
									_instance.B0[d].V_values.push_back(sign * kernels[master].dense_values[kernels[master].rows * c + _DOFs * (masterit->DOFOffset + nodes[n] - masterit->begin) + dof]);
								}
							}
						}
					} else {
						auto dit = _mesh.nodes->dintervals[d].begin();
						for (size_t n = _mesh.FETIData->inodesDistribution[i]; n < _mesh.FETIData->inodesDistribution[i + 1]; n++) {
							while (dit->end < nodes[n]) {
								++dit;
							}
							for (int dof = 0; dof < _DOFs; dof++) {
								_instance.B0[d].I_row_indices.push_back(rowIndex[i] + 1);
								_instance.B0[d].J_col_indices.push_back(_DOFs * (dit->DOFOffset + nodes[n] - dit->begin) + dof + 1);
								_instance.B0[d].V_values.push_back(sign);
							}
						}
					}
				}
			}
			_instance.B0[d].rows = rCounters[_mesh.elements->clusters[d]];
			_instance.B0[d].cols = _instance.domainDOFCount[d];
			_instance.B0[d].nnz = _instance.B0[d].I_row_indices.size();
			_instance.B0subdomainsMap[d].reserve(_instance.B0[d].nnz);
			for (eslocal i = _instance.B0subdomainsMap[d].size(); i < _instance.B0[d].nnz; i++) {
				_instance.B0subdomainsMap[d].push_back(_instance.B0[d].I_row_indices[i]);
			}
		}
	}
}

void EqualityConstraints::B0Corners()
{
	for (size_t d = 0; d < _instance.domains; d++) {
		_instance.B0[d].cols = _instance.K[d].cols;
	}

	size_t lambdas = 1;

	auto domains = _mesh.FETIData->cornerDomains->cbegin();
	for (size_t n = 0; n < _mesh.FETIData->corners.size(); ++n, ++domains) {
		for (size_t dof = 0; dof < 1; dof++) {
			for (size_t d1 = 0, d2 = 1; d2 < domains->size(); ++d1, ++d2) {

				auto d1it = _mesh.nodes->dintervals[domains->at(d1)].begin();
				auto d2it = _mesh.nodes->dintervals[domains->at(d2)].begin();
				while (d1it->end < _mesh.FETIData->corners[n]) {
					++d1it;
				}
				while (d2it->end < _mesh.FETIData->corners[n]) {
					++d2it;
				}

				for (int dof = 0; dof < _DOFs; dof++) {
					_instance.B0[domains->at(d1)].I_row_indices.push_back(lambdas);
					_instance.B0[domains->at(d1)].J_col_indices.push_back(_DOFs * (_mesh.FETIData->corners[n] - d1it->begin + d1it->DOFOffset) + dof + 1);
					_instance.B0[domains->at(d1)].V_values.push_back(1);

					_instance.B0[domains->at(d2)].I_row_indices.push_back(lambdas);
					_instance.B0[domains->at(d2)].J_col_indices.push_back(_DOFs * (_mesh.FETIData->corners[n] - d2it->begin + d2it->DOFOffset) + dof + 1);
					_instance.B0[domains->at(d2)].V_values.push_back(-1);

					lambdas++;
				}
			}
		}
	}

	#pragma omp parallel for
	for  (size_t p = 0; p < _instance.domains; p++) {
		_instance.B0[p].rows = lambdas - 1;
		_instance.B0[p].cols = _instance.domainDOFCount[p];
		_instance.B0[p].nnz = _instance.B0[p].I_row_indices.size();

		_instance.B0subdomainsMap[p].reserve(_instance.B0[p].nnz);
		for (eslocal i = _instance.B0subdomainsMap[p].size(); i < _instance.B0[p].nnz; i++) {
			_instance.B0subdomainsMap[p].push_back(_instance.B0[p].I_row_indices[i] - 1);
		}
	}
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
	std::vector<int> interfaceNodes;

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

