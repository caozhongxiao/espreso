
#include "constraints.h"

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
#include "../../config/ecf/physics/physics.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

Constraints::Constraints(Instance &instance, Mesh &mesh, const RegionMap<ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling)
: _instance(instance), _mesh(mesh), _contactSize(0), _DOFs(1)
{
	update(dirichlet, withRedundantMultiplier, withScaling);
}

Constraints::Constraints(Instance &instance, Mesh &mesh, const RegionMap<ECFExpressionOptionalVector> &dirichlet, int DOFs, bool withRedundantMultiplier, bool withScaling)
: _instance(instance), _mesh(mesh), _contactSize(0), _DOFs(DOFs)
{
	update(dirichlet, withRedundantMultiplier, withScaling);
}

eslocal Constraints::computeIntervalsOffsets(std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, eslocal)> setsize)
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

void Constraints::update(const RegionMap<ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling)
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

void Constraints::update(const RegionMap<ECFExpressionOptionalVector> &dirichlet, bool withRedundantMultiplier, bool withScaling)
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

void Constraints::update(bool withRedundantMultiplier, bool withScaling)
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

