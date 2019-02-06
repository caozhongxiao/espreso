
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "physics/assembler/dataholder.h"
#include "uniformnodesfeticomposer.h"

#include "physics/assembler/controllers/controller.h"
#include "physics/assembler/provider/feti/fetiprovider.h"

#include "physics/assembler/assembler.h"
#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"
#include "basis/logging/logging.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"

#include "config/ecf/linearsolver/feti.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/surfacestore.h"

#include "solver/generic/SparseMatrix.h"

#include "wrappers/bem/bemwrapper.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void UniformNodesFETIComposer::initDOFs()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	_domainDOFsSize.resize(info::mesh->elements->ndomains);

	// nID, domain
	std::vector<std::vector<std::pair<esint, esint> > > ntodomains(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<esint, esint> > tdata;

		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			std::vector<esint> dnodes(
					(info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d])->begin(),
					(info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d + 1])->begin());

			Esutils::sortAndRemoveDuplicity(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<esint, esint>(dnodes[i], d));
			}
			_domainDOFsSize[d] = dnodes.size() * _DOFs;
		}

		ntodomains[t].swap(tdata);
	}

	for (size_t t = 1; t < threads; t++) {
		ntodomains[0].insert(ntodomains[0].end(), ntodomains[t].begin(), ntodomains[t].end());
	}

	std::sort(ntodomains[0].begin(), ntodomains[0].end());

	std::vector<std::vector<esint> > DOFs(threads);
	std::vector<size_t> ndistribution = tarray<size_t>::distribute(threads, ntodomains[0].size());

	for (size_t t = 1; t < threads; t++) {
		while (
				ndistribution[t] < ntodomains[0].size() &&
				ntodomains[0][ndistribution[t]].first == ntodomains[0][ndistribution[t] - 1].first) {

			++ndistribution[t];
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tDOFs(info::mesh->elements->ndomains);

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			tDOFs[ntodomains[0][n].second] += _DOFs;
		}

		DOFs[t].swap(tDOFs);
	}

	Esutils::sizesToOffsets(DOFs);

	std::vector<std::vector<std::vector<esint> > > sBuffer(threads);
	std::vector<std::vector<esint> > rBuffer(info::mesh->neighboursWithMe.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nranks = info::mesh->nodes->ranks->begin() + ntodomains[0][ndistribution[t]].first;

		std::vector<std::vector<esint> > tBuffer(info::mesh->neighboursWithMe.size());

		size_t n = ndistribution[t];
		while (n < ndistribution[t + 1]) {
			size_t begin = n++;
			while (n < ntodomains[0].size() && ntodomains[0][n].first == ntodomains[0][n - 1].first) {
				++n;
			}

			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (info::mesh->neighboursWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(n - begin);
				for (size_t i = begin; i < n; i++) {
					tBuffer[noffset].push_back(info::mesh->elements->firstDomain + ntodomains[0][i].second);
					for (int dof = 0; dof < _DOFs; ++dof) {
						tBuffer[noffset].push_back(DOFs[t][ntodomains[0][i].second] + dof);
					}
				}
			}
			++nranks;

			for (size_t i = begin; i < n; i++) {
				DOFs[t][ntodomains[0][i].second] += _DOFs;
			}
		}
		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < sBuffer[t].size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, info::mesh->neighboursWithMe)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	// domain0, DOF0, DOF1, ..., DOFn, domain1, DOF0, DOF1, ..., DOFn, ...; domain0, ...
	std::vector<esint> DOFDistribution(1), DOFData;

	// TODO: make it parallel
	// parallelization is possible if node order will be kept as: boundary first!
	// now we prefer generality
	auto nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> roffset(rBuffer.size());
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		esint noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (info::mesh->neighboursWithMe[noffset] < *r) {
				++noffset;
			}

			esint domains = rBuffer[noffset][roffset[noffset]++];
			for (esint d = 0; d < domains; ++d) {
				DOFData.push_back(rBuffer[noffset][roffset[noffset]++]);
				for (int dof = 0; dof < _DOFs; ++dof) {
					DOFData.push_back(rBuffer[noffset][roffset[noffset]++]);
				}
			}
		}
		DOFDistribution.push_back(DOFData.size());
	}

	std::vector<size_t> distribution = info::mesh->nodes->distribution, datadistribution(threads + 1);
	for (size_t t = 1; t < threads; t++) {
		++distribution[t];
		datadistribution[t] = DOFDistribution[distribution[t]];
	}
	datadistribution[threads] = DOFDistribution[distribution[threads]];
	++distribution[threads];


	_DOFMap = new serializededata<esint, esint>(
			tarray<esint>(distribution, DOFDistribution),
			tarray<esint>(datadistribution, DOFData));
}

void UniformNodesFETIComposer::buildDirichlet()
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

void UniformNodesFETIComposer::buildPatterns()
{
	_KPermutation.resize(info::mesh->elements->ndomains);
	_RHSPermutation.resize(info::mesh->elements->ndomains);

	data->K.resize(info::mesh->elements->ndomains);
	data->M.resize(info::mesh->elements->ndomains);
	data->f.resize(info::mesh->elements->ndomains);
	data->R.resize(info::mesh->elements->ndomains);
	data->primalSolution.resize(info::mesh->elements->ndomains);


	size_t threads = info::env::OMP_NUM_THREADS;
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			if (_BEMDomain[d]) {
				buildKBEMPattern(d);
			} else {
				buildKFEMPattern(d);
			}
		}
	}

	buildB1Pattern();

	_foreignDOFs = _DOFs * (info::mesh->nodes->size - info::mesh->nodes->uniqueSize);
}

void UniformNodesFETIComposer::buildKFEMPattern(esint domain)
{
	data->K[domain].rows = _domainDOFsSize[domain];
	data->K[domain].cols = _domainDOFsSize[domain];
	data->f[domain].resize(_domainDOFsSize[domain]);
	data->R[domain].resize(_domainDOFsSize[domain]);

	MatrixType mtype = _provider.getMatrixType(domain);

	auto ebegin = info::mesh->elements->procNodes->cbegin() + info::mesh->elements->elementsDistribution[domain];
	auto eend = info::mesh->elements->procNodes->cbegin() + info::mesh->elements->elementsDistribution[domain + 1];

	esint Ksize = 0, RHSsize = 0;
	for (auto e = ebegin; e != eend; ++e) {
		RHSsize += e->size() * _DOFs;
		Ksize += getMatrixSize(e->size() * _DOFs, mtype);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (
				info::mesh->boundaryRegions[r]->dimension &&
				info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {

			esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
			esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
			auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

			for (esint i = begin; i < end; ++i, ++enodes) {
				RHSsize += enodes->size() * _DOFs;
				Ksize += getMatrixSize(enodes->size() * _DOFs, mtype);;
			}
		}
	}

	std::vector<esint> permK(Ksize), permRHS(RHSsize);
	std::vector<IJ> KPattern(Ksize);
	std::vector<esint> RHSPattern(RHSsize), ROW, COL;

	IJ *Koffset = KPattern.data();
	esint *RHSoffset = RHSPattern.data();

	auto insert = [&] (serializededata<esint, esint>::const_iterator &enodes) {
		esint *_RHS = RHSoffset;
		for (auto n = enodes->begin(); n != enodes->end(); ++n, ++RHSoffset) {
			auto DOFs = (_DOFMap->begin() + *n)->begin();
			while (*DOFs != domain + info::mesh->elements->firstDomain) {
				DOFs += 1 + _DOFs;
			}
			*RHSoffset = *(DOFs + 1);
		}
		for (int dof = 1; dof < _DOFs; ++dof) {
			for (size_t n = 0; n < enodes->size(); ++n, ++RHSoffset) {
				*RHSoffset = *(_RHS + n) + dof;
			}
		}
		insertKPattern(Koffset, _RHS, RHSoffset, mtype);
	};

	for (auto e = ebegin; e != eend; ++e) {
		insert(e);
		Koffset += getMatrixSize(e->size() * _DOFs, mtype);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (
				info::mesh->boundaryRegions[r]->dimension &&
				info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {

			esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
			esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
			auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

			for (esint i = begin; i < end; ++i, ++enodes) {
				insert(enodes);
				Koffset += getMatrixSize(enodes->size() * _DOFs, mtype);
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

	ROW.reserve(_domainDOFsSize[domain] + 1);
	COL.reserve(KPattern.size());
	ROW.push_back(1);
	COL.push_back(KPattern[pK.front()].column + 1);
	permK[pK.front()] = 0;
	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
			++nonzeros;
			COL.push_back(KPattern[pK[i]].column + 1);
			if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
				ROW.push_back(nonzeros + 1);
			}
		}
		permK[pK[i]] = nonzeros;
	}
	permRHS[pRHS.front()] = 0;
	for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
			++nonzeros;
		}
		permRHS[pRHS[i]] = nonzeros;
	}

	_KPermutation[domain].swap(permK);
	_RHSPermutation[domain].swap(permRHS);

	ROW.push_back(COL.size() + 1);

	data->K[domain].nnz = COL.size();
	data->K[domain].mtype = _provider.getMatrixType(domain);
	switch (data->K[domain].mtype) {
	case MatrixType::REAL_UNSYMMETRIC: data->K[domain].type = 'G'; break;
	default: data->K[domain].type = 'S';
	}
	data->K[domain].CSR_V_values.resize(COL.size());
	data->K[domain].CSR_I_row_indices.swap(ROW);
	data->K[domain].CSR_J_col_indices.swap(COL);
	data->M[domain] = data->K[domain];
	data->M[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	data->M[domain].type = 'S';
}
void UniformNodesFETIComposer::buildKBEMPattern(esint domain)
{
	SurfaceStore* surface = info::mesh->domainsSurface;
	data->K[domain].rows = surface->cdistribution[domain + 1] - surface->cdistribution[domain];
	data->K[domain].cols = surface->cdistribution[domain + 1] - surface->cdistribution[domain];
	data->K[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	data->K[domain].type = 'S';
	data->f[domain].resize(data->K[domain].rows);

	data->K[domain].CSR_I_row_indices.push_back(1);
	for (esint r = 0; r < data->K[domain].rows; r++) {
		for (esint c = r; c < data->K[domain].cols; c++) {
			data->K[domain].CSR_J_col_indices.push_back(c + 1);
		}
		data->K[domain].CSR_I_row_indices.push_back(data->K[domain].CSR_J_col_indices.size() + 1);
	}

	data->K[domain].CSR_V_values.resize(data->K[domain].CSR_J_col_indices.size());
	data->K[domain].nnz = data->K[domain].CSR_J_col_indices.size();

	esint RHSsize = 0;
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (
				info::mesh->boundaryRegions[r]->dimension &&
				info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {

			esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
			esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
			auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

			for (esint i = begin; i < end; ++i, ++enodes) {
				RHSsize += enodes->size() * _DOFs;
			}
		}
	}

	std::vector<esint> permRHS(RHSsize);
	std::vector<esint> RHSPattern(RHSsize), ROW, COL;

	esint *RHSoffset = RHSPattern.data();

	auto insert = [&] (serializededata<esint, esint>::const_iterator &enodes) {
		esint *_RHS = RHSoffset;
		for (auto n = enodes->begin(); n != enodes->end(); ++n, ++RHSoffset) {
			auto DOFs = (_DOFMap->begin() + *n)->begin();
			while (*DOFs != domain + info::mesh->elements->firstDomain) {
				DOFs += 1 + _DOFs;
			}
			*RHSoffset = *(DOFs + 1);
		}
		for (int dof = 1; dof < _DOFs; ++dof) {
			for (size_t n = 0; n < enodes->size(); ++n, ++RHSoffset) {
				*RHSoffset = *(_RHS + n) + dof;
			}
		}
	};

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (
				info::mesh->boundaryRegions[r]->dimension &&
				info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {

			esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
			esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
			auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

			for (esint i = begin; i < end; ++i, ++enodes) {
				insert(enodes);
			}
		}
	}

	std::vector<esint> pRHS(RHSPattern.size());
	std::iota(pRHS.begin(), pRHS.end(), 0);
	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	if (pRHS.size()) {
		permRHS[pRHS.front()] = 0;
		for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
			if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
				++nonzeros;
			}
			permRHS[pRHS[i]] = nonzeros;
		}
	}

	_RHSPermutation[domain].swap(permRHS);
}

void UniformNodesFETIComposer::buildB1Pattern()
{
	esint doffset = 0;

	data->B1.resize(info::mesh->elements->ndomains);
	data->B1c.resize(info::mesh->elements->ndomains);
	data->B1subdomainsMap.resize(info::mesh->elements->ndomains);
	data->B1duplicity.resize(info::mesh->elements->ndomains);
	data->LB.resize(info::mesh->elements->ndomains);
	data->block.resize(3);

	auto dmap = _DOFMap->begin();
	for (size_t i = 0, prev = 0; i < _dirichletMap.size(); prev = _dirichletMap[i++] / _DOFs) {
		while (i + 1 < _dirichletMap.size() && _dirichletMap[i + 1] == _dirichletMap[i]) { ++i; }
		dmap += (_dirichletMap[i] / _DOFs) - prev;
		for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				++doffset;
			}
			if (!_configuration.redundant_lagrange) {
				break;
			}
		}
	}

	esint dsize = Communication::exscan(doffset);

	dmap = _DOFMap->begin();
	for (size_t i = 0, prev = 0; i < _dirichletMap.size(); prev = _dirichletMap[i++] / _DOFs) {
		while (i + 1 < _dirichletMap.size() && _dirichletMap[i + 1] == _dirichletMap[i]) { ++i; }
		dmap += (_dirichletMap[i] / _DOFs) - prev;
		for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				data->B1[*d - info::mesh->elements->firstDomain].I_row_indices.push_back(doffset + 1);
				data->B1[*d - info::mesh->elements->firstDomain].J_col_indices.push_back(*(d + _dirichletMap[i] % _DOFs + 1) + 1);
				data->B1clustersMap.push_back({ doffset, info::mpi::rank });
				++doffset;
			}
			if (!_configuration.redundant_lagrange) {
				break;
			}
		}
	}

	_domainDirichletSize.resize(info::mesh->elements->ndomains);
	for (esint d = 0; d < info::mesh->elements->ndomains; ++d) {
		data->B1[d].V_values.resize(data->B1[d].I_row_indices.size(), 1);
		data->B1c[d].resize(data->B1[d].I_row_indices.size());
		data->B1duplicity[d].resize(data->B1[d].I_row_indices.size(), 1);
		_domainDirichletSize[d] = data->B1[d].I_row_indices.size();
	}

	data->block[DataHolder::CONSTRAINT::DIRICHLET] = dsize;

	esint goffset = 0;

	dmap = _DOFMap->begin();
	auto nranks = info::mesh->nodes->ranks->begin();
	auto exclude = _dirichletMap.begin();
	std::vector<std::vector<esint> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	auto send = [&] () {
		esint noffset = 0;
		for (auto r = nranks->begin() + 1; r != nranks->end(); ++r) {
			while (info::mesh->neighbours[noffset] < *r) {
				++noffset;
			}
			sBuffer[noffset].push_back(goffset);
		}
	};

	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (dmap->size() / (1 + _DOFs) > 1) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				esint ndomains = dmap->size() / (1 + _DOFs);
				if (_configuration.redundant_lagrange) {
					while (exclude != _dirichletMap.end() && *exclude < n * _DOFs + dof) {
						++exclude;
					}
					if (exclude == _dirichletMap.end() || n * _DOFs + dof != *exclude) {
						if (*nranks->begin() == info::mpi::rank) {
							send();
							goffset += ndomains * (ndomains - 1) / 2;
						}
					}
				} else {
					if (*nranks->begin() == info::mpi::rank) {
						send();
						goffset += ndomains - 1;
					}
				}
			}
		}
	}

	esint gsize = Communication::exscan(goffset);

	for (size_t n = 0; n < info::mesh->neighbours.size(); ++n) {
		for (size_t i = 0; i < sBuffer[n].size(); ++i) {
			sBuffer[n][i] += goffset;
		}
	}

	if (!Communication::receiveLowerUnknownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange gluing offsets";
	}

	std::vector<esint> dDistribution = info::mesh->elements->gatherDomainsProcDistribution();

	dmap = _DOFMap->begin();
	nranks = info::mesh->nodes->ranks->begin();
	exclude = _dirichletMap.begin();

	auto push = [&] (esint d, esint lambda, esint dof, double value, esint domains) {
		d -= info::mesh->elements->firstDomain;
		data->B1[d].I_row_indices.push_back(lambda + 1);
		data->B1[d].J_col_indices.push_back(dof + 1);
		data->B1[d].V_values.push_back(value);
		data->B1duplicity[d].push_back(1. / domains);
	};

	auto fill = [&] (esint d1, esint d2, esint lambda, esint dof1, esint dof2, esint domains) {
		if (info::mesh->elements->firstDomain <= d1 && d1 < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
			push(d1, lambda, dof1, 1, domains);
		}
		if (info::mesh->elements->firstDomain <= d2 && d2 < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
			push(d2, lambda, dof2, -1, domains);
		}
		if (
				(info::mesh->elements->firstDomain <= d1 && d1 < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) ||
				(info::mesh->elements->firstDomain <= d2 && d2 < info::mesh->elements->firstDomain + info::mesh->elements->ndomains)) {

			data->B1clustersMap.push_back({ lambda, info::mpi::rank });

			if (d1 <info::mesh->elements->firstDomain || info::mesh->elements->firstDomain + info::mesh->elements->ndomains <= d1) {
				data->B1clustersMap.back().push_back(std::lower_bound(dDistribution.begin(), dDistribution.end(), d1 + 1) - dDistribution.begin() - 1);
			}
			if (d2 <info::mesh->elements->firstDomain || info::mesh->elements->firstDomain + info::mesh->elements->ndomains <= d2) {
				data->B1clustersMap.back().push_back(std::lower_bound(dDistribution.begin(), dDistribution.end(), d2 + 1) - dDistribution.begin() - 1);
			}
		}
	};

	std::vector<esint> roffset(info::mesh->neighbours.size());
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (dmap->size() / (1 + _DOFs) > 1) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				if (_configuration.redundant_lagrange) {
					while (exclude != _dirichletMap.end() && *exclude < n * _DOFs + dof) {
						++exclude;
					}
				}
				esint lambda = goffset;
				esint ndomains = dmap->size() / (1 + _DOFs);
				if (*nranks->begin() != info::mpi::rank) {
					if (!_configuration.redundant_lagrange || exclude == _dirichletMap.end() || n * _DOFs + dof != *exclude) {
						esint noffset = 0;
						while (info::mesh->neighbours[noffset] < *nranks->begin()) {
							++noffset;
						}
						lambda = rBuffer[noffset][roffset[noffset]++];
					}
				}
				if (_configuration.redundant_lagrange) {
					if (exclude == _dirichletMap.end() || n * _DOFs + dof != *exclude) {
						for (auto d1 = dmap->begin(); d1 != dmap->end(); d1 += 1 + _DOFs) {
							for (auto d2 = d1 + 1 + _DOFs; d2 != dmap->end(); d2 += 1 + _DOFs, ++lambda) {
								fill(*d1, *d2, dsize + lambda, *(d1 + 1 + dof), *(d2 + 1 + dof), ndomains);
							}
						}
						if (*nranks->begin() == info::mpi::rank) {
							goffset += ndomains * (ndomains - 1) / 2;
						}
					}
				} else {
					for (auto d1 = dmap->begin(), d2 = dmap->begin() + 1 + _DOFs; d2 != dmap->end(); d1 = d2, d2 += 1 + _DOFs, ++lambda) {
						fill(*d1, *d2, dsize + lambda, *(d1 + 1 + dof), *(d2 + 1 + dof), ndomains);
					}
					if (*nranks->begin() == info::mpi::rank) {
						goffset += ndomains - 1;
					}
				}
			}
		}
	}

	for (esint d = 0; d < info::mesh->elements->ndomains; ++d) {
		data->B1c[d].resize(data->B1[d].I_row_indices.size());

		data->B1[d].nnz = data->B1[d].I_row_indices.size();
		data->B1[d].rows = dsize + gsize;
		data->B1[d].cols = data->K[d].rows;
		data->LB[d].resize(data->B1[d].nnz, -std::numeric_limits<double>::infinity());

		data->B1subdomainsMap[d] = data->B1[d].I_row_indices;
		for (size_t i = 0; i < data->B1subdomainsMap[d].size(); ++i) {
			--data->B1subdomainsMap[d][i];
		}
	}

	data->block[DataHolder::CONSTRAINT::EQUALITY_CONSTRAINTS] = dsize + gsize;
}


void UniformNodesFETIComposer::assemble(Matrices matrices, const SolverParameters &parameters)
{
	if (!(matrices & (Matrices::K | Matrices::M | Matrices::R | Matrices::f))) {
		return;
	}

	#pragma omp parallel for
	for  (esint d = 0; d < info::mesh->elements->ndomains; d++) {

		size_t KIndex = 0, RHSIndex = 0;
		double KReduction = parameters.timeIntegrationConstantK, RHSReduction = parameters.internalForceReduction;
		Controller::InstanceFiller filler;

		switch (_provider.getMatrixType(d)) {
		case MatrixType::REAL_UNSYMMETRIC:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						data->f[d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						data->R[d][_RHSPermutation[d][RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = 0; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							data->K[d].CSR_V_values[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							data->M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						data->f[d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						data->R[d][_RHSPermutation[d][RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							data->K[d].CSR_V_values[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							data->M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		}

		clearMatrices(matrices, d);
		filler.begin = info::mesh->elements->elementsDistribution[d];
		filler.end = info::mesh->elements->elementsDistribution[d + 1];

		if (_BEMDomain[d]) {
			_controler.processBEMdomain(d, data->K[d].CSR_V_values.data());
		} else {
			_controler.processElements(matrices, parameters, filler);
		}

		KReduction = parameters.internalForceReduction;
		filler.Me.resize(0, 0);
		filler.Re.resize(0, 0);

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->distribution.size()) {
				if (info::mesh->boundaryRegions[r]->eintervalsDistribution[d] < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]) {
					filler.begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[d]].begin;
					filler.end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					_controler.processBoundary(matrices, parameters, r, filler);
				}
			}
		}
	};
}

void UniformNodesFETIComposer::setDirichlet(Matrices matrices, double reduction, const std::vector<double> &subtraction)
{
	std::vector<double> values(_dirichletMap.size());
	_controler.dirichletValues(values);

	for (size_t i = 0; i < values.size(); i++) {
		values[i] *= reduction;
	}

	if (subtraction.size()) {
		for (size_t i = 0; i < _dirichletMap.size(); i++) {
			values[_dirichletPermutation[i]] -= subtraction[_dirichletMap[i]];
		}
	}

	std::vector<esint> doffset(info::mesh->elements->ndomains);

	auto dmap = _DOFMap->begin();
	for (size_t i = 0, j = 0, prev = 0; i < _dirichletMap.size(); prev = _dirichletMap[i] / _DOFs, i = j) {
		double value = 0;
		while (j < _dirichletMap.size() && _dirichletMap[j] == _dirichletMap[i]) {
			value += values[_dirichletPermutation[j++]];
		}
		value /= j - i;
		dmap += (_dirichletMap[i] / _DOFs) - prev;
		for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				data->B1c[*d - info::mesh->elements->firstDomain][doffset[*d - info::mesh->elements->firstDomain]++] = value;
			}
			if (!_configuration.redundant_lagrange) {
				break;
			}
		}
	}
	if (_configuration.scaling && _configuration.redundant_lagrange) {
		updateDuplicity();
	}
}


void UniformNodesFETIComposer::updateDuplicity()
{
	std::vector<std::vector<double> > diagonals(info::mesh->elements->ndomains);
	#pragma omp parallel for
	for  (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		diagonals[d] = data->K[d].getDiagonal();
	}

	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	auto dmap = _DOFMap->begin();
	auto nranks = info::mesh->nodes->ranks->begin();
	std::vector<double> buffer;

	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (nranks->size() > 1) {
			buffer.clear();
			for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
				if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
					for (int dof = 0; dof < _DOFs; ++dof) {
						buffer.push_back(diagonals[*d - info::mesh->elements->firstDomain][*(d + 1 + dof)]);
					}
				}
			}

			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbours[noffset] < *r) {
						++noffset;
					}
					sBuffer[noffset].insert(sBuffer[noffset].end(), buffer.begin(), buffer.end());
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange diagonal values";
	}

	dmap = _DOFMap->begin();
	nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> roffset(info::mesh->neighbours.size());
	std::vector<esint> dDistribution = info::mesh->elements->gatherDomainsProcDistribution();
	auto exclude = _dirichletMap.begin();

	auto push = [&] (esint d1, esint d2, double v1, double v2, double sum) {
		if (info::mesh->elements->firstDomain <= d1 && d1 < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
			data->B1duplicity[d1 - info::mesh->elements->firstDomain].push_back(v2 / sum);
		}
		if (info::mesh->elements->firstDomain <= d2 && d2 < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
			data->B1duplicity[d2 - info::mesh->elements->firstDomain].push_back(v1 / sum);
		}
	};

	for (size_t d = 0; d < data->B1duplicity.size(); d++) {
		data->B1duplicity[d].resize(_domainDirichletSize[d]);
	}
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (dmap->size() / (1 + _DOFs) > 1) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				buffer.clear();

				for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
					if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
						buffer.push_back(diagonals[*d - info::mesh->elements->firstDomain][*(d + 1 + dof)]);
					} else {
						esint r = std::lower_bound(dDistribution.begin(), dDistribution.end(), *d + 1) - dDistribution.begin() - 1;
						esint noffset = 0;
						while (info::mesh->neighbours[noffset] < r) {
							++noffset;
						}
						buffer.push_back(rBuffer[noffset][roffset[noffset]++]);
					}
				}

				double sum = 0;
				for (size_t i = 0; i < buffer.size(); i++) {
					sum += buffer[i];
				}

				while (exclude != _dirichletMap.end() && *exclude < n * _DOFs + dof) {
					++exclude;
				}
				if (exclude == _dirichletMap.end() || n * _DOFs + dof != *exclude) {
					esint v1 = 0, v2 = 1;
					for (auto d1 = dmap->begin(); d1 != dmap->end(); d1 += 1 + _DOFs, ++v1, v2 = v1 + 1) {
						for (auto d2 = d1 + 1 + _DOFs; d2 != dmap->end(); d2 += 1 + _DOFs, ++v2) {
							push(*d1, *d2, buffer[v1], buffer[v2], sum);
						}
					}
				}
			}
		}
	}
}

void UniformNodesFETIComposer::fillSolution()
{
	for (size_t i = 0; i < _BEMDomain.size(); i++) {
		if (_BEMDomain[i]) {
			data->primalSolution[i].resize(_domainDOFsSize[i]);
			_controler.fillBEMInterior(i, data->primalSolution[i].data());
		}
	}
	avgGather(_controler.solution()->data, data->primalSolution);
}

void UniformNodesFETIComposer::_divide(std::vector<double> &in, std::vector<std::vector<double> > &out, bool split)
{
	out.resize(info::mesh->elements->ndomains);
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		out[d].resize(_domainDOFsSize[d]);
	}
	size_t i = 0;
	esint dbegin = info::mesh->elements->firstDomain;
	esint dend = info::mesh->elements->firstDomain + info::mesh->elements->ndomains;
	for (auto n = _DOFMap->begin(); n != _DOFMap->end(); ++n, ++i) {
		esint domains = n->size() / (_DOFs + 1);
		for (esint d = 0; d < domains; ++d) {
			for (esint dof = 0; dof < _DOFs; ++dof) {
				esint domain = n->at(d * (_DOFs + 1));
				esint index = n->at(d * (_DOFs + 1) + dof + 1);
				if (dbegin <= domain && domain < dend) {
					if (split) {
						out[domain - dbegin][index] = in[i] / domains;
					} else {
						out[domain - dbegin][index] = in[i];
					}
				}
			}
		}
	}
}

void UniformNodesFETIComposer::_gather(std::vector<double> &out, std::vector<std::vector<double> > &in, bool avg)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	out.resize(_DOFs * info::mesh->nodes->size);

	std::vector<std::vector<std::vector<double> > > sBuffer(threads);
	std::vector<std::vector<double> > rBuffer(info::mesh->neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t i = info::mesh->nodes->distribution[t];
		auto nranks = info::mesh->nodes->ranks->begin(t);

		std::vector<std::vector<double> > tBuffer(info::mesh->neighbours.size());

		for (auto map = _DOFMap->begin(t); map != _DOFMap->end(t); ++map, ++i, ++nranks) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				out[i * _DOFs + dof] = 0;
			}
			for (auto d = map->begin(); d != map->end(); d += 1 + _DOFs) {
				if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
					for (int dof = 0; dof < _DOFs; ++dof) {
						double value = in[*d - info::mesh->elements->firstDomain][*(d + 1 + dof)];
						if (avg) {
							value /= ((map->end() - map->begin()) / (1 + _DOFs));
						}
						out[i * _DOFs + dof] += value;
					}
				}
			}

			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbours[noffset] < *r) {
						++noffset;
					}

					for (esint dof = 0; dof < _DOFs; ++dof) {
						tBuffer[noffset].push_back(out[i * _DOFs + dof]);
					}
				}
			}
		}

		sBuffer[t].swap(tBuffer);
	}

	for (size_t n = 0; n < sBuffer[0].size(); ++n) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
		rBuffer[n].resize(sBuffer[0][n].size());
	}

	if (!Communication::exchangeKnownSize(sBuffer[0], rBuffer, info::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize solution.";
	}

	std::vector<esint> roffset(info::mesh->neighbours.size());
	auto nranks = info::mesh->nodes->ranks->begin();
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		esint noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			if (*r != info::mpi::rank) {
				while (info::mesh->neighbours[noffset] < *r) {
					++noffset;
				}

				for (esint dof = 0; dof < _DOFs; ++dof) {
					out[n * _DOFs + dof] += rBuffer[noffset][roffset[noffset]++];
				}
			}
		}
	}
}

