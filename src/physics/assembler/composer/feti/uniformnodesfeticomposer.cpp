
#include "uniformnodesfeticomposer.h"

#include "../../controllers/controller.h"
#include "../../provider/feti/fetiprovider.h"

#include "../../../dataholder.h"

#include "../../../../globals/run.h"
#include "../../../../basis/containers/serializededata.h"
#include "../../../../basis/matrices/denseMatrix.h"
#include "../../../../basis/logging/logging.h"
#include "../../../../basis/utilities/utils.h"
#include "../../../../basis/utilities/communication.h"

#include "../../../../config/ecf/environment.h"
#include "../../../../config/ecf/solver/feti.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"
#include "../../../../mesh/store/nodestore.h"

#include "../../../../solver/generic/SparseMatrix.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void UniformNodesFETIComposer::initDOFs()
{
	size_t threads = environment->OMP_NUM_THREADS;

	_domainDOFsSize.resize(run::mesh->elements->ndomains);

	// nID, domain
	std::vector<std::vector<std::pair<eslocal, eslocal> > > ntodomains(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<eslocal, eslocal> > tdata;

		for (eslocal d = run::mesh->elements->domainDistribution[t]; d != run::mesh->elements->domainDistribution[t + 1]; ++d) {
			std::vector<eslocal> dnodes(
					(run::mesh->elements->procNodes->begin() + run::mesh->elements->elementsDistribution[d])->begin(),
					(run::mesh->elements->procNodes->begin() + run::mesh->elements->elementsDistribution[d + 1])->begin());

			Esutils::sortAndRemoveDuplicity(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<eslocal, eslocal>(dnodes[i], d));
			}
			_domainDOFsSize[d] = dnodes.size() * _DOFs;
		}

		ntodomains[t].swap(tdata);
	}

	for (size_t t = 1; t < threads; t++) {
		ntodomains[0].insert(ntodomains[0].end(), ntodomains[t].begin(), ntodomains[t].end());
	}

	std::sort(ntodomains[0].begin(), ntodomains[0].end());

	std::vector<std::vector<eslocal> > DOFs(threads);
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
		std::vector<eslocal> tDOFs(run::mesh->elements->ndomains);

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			tDOFs[ntodomains[0][n].second] += _DOFs;
		}

		DOFs[t].swap(tDOFs);
	}

	Esutils::sizesToOffsets(DOFs);

	std::vector<std::vector<std::vector<eslocal> > > sBuffer(threads);
	std::vector<std::vector<eslocal> > rBuffer(run::mesh->neighboursWithMe.size());


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nranks = run::mesh->nodes->ranks->begin() + ntodomains[0][ndistribution[t]].first;

		std::vector<std::vector<eslocal> > tBuffer(run::mesh->neighboursWithMe.size());

		size_t n = ndistribution[t];
		while (n < ndistribution[t + 1]) {
			size_t begin = n++;
			while (n < ntodomains[0].size() && ntodomains[0][n].first == ntodomains[0][n - 1].first) {
				++n;
			}

			eslocal noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (run::mesh->neighboursWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(n - begin);
				for (size_t i = begin; i < n; i++) {
					tBuffer[noffset].push_back(run::mesh->elements->firstDomain + ntodomains[0][i].second);
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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, run::mesh->neighboursWithMe)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	// domain0, DOF0, DOF1, ..., DOFn, domain1, DOF0, DOF1, ..., DOFn, ...; domain0, ...
	std::vector<eslocal> DOFDistribution(1), DOFData;

	// TODO: make it parallel
	// parallelization is possible if node order will be kept as: boundary first!
	// now we prefer generality
	auto nranks = run::mesh->nodes->ranks->begin();
	std::vector<eslocal> roffset(rBuffer.size());
	for (eslocal n = 0; n < run::mesh->nodes->size; ++n, ++nranks) {
		eslocal noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (run::mesh->neighboursWithMe[noffset] < *r) {
				++noffset;
			}

			eslocal domains = rBuffer[noffset][roffset[noffset]++];
			for (eslocal d = 0; d < domains; ++d) {
				DOFData.push_back(rBuffer[noffset][roffset[noffset]++]);
				for (int dof = 0; dof < _DOFs; ++dof) {
					DOFData.push_back(rBuffer[noffset][roffset[noffset]++]);
				}
			}
		}
		DOFDistribution.push_back(DOFData.size());
	}

	std::vector<size_t> distribution = run::mesh->nodes->distribution, datadistribution(threads + 1);
	for (size_t t = 1; t < threads; t++) {
		++distribution[t];
		datadistribution[t] = DOFDistribution[distribution[t]];
	}
	datadistribution[threads] = DOFDistribution[distribution[threads]];
	++distribution[threads];


	_DOFMap = new serializededata<eslocal, eslocal>(
			tarray<eslocal>(distribution, DOFDistribution),
			tarray<eslocal>(datadistribution, DOFData));
}

void UniformNodesFETIComposer::initDirichlet()
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

void UniformNodesFETIComposer::buildPatterns()
{
	buildKPattern();
	buildB1Pattern();
	buildB0Pattern();
}

void UniformNodesFETIComposer::buildKPattern()
{
	size_t threads = environment->OMP_NUM_THREADS;

	_KPermutation.resize(run::mesh->elements->ndomains);
	_RHSPermutation.resize(run::mesh->elements->ndomains);

	run::data->K.resize(run::mesh->elements->ndomains);
	run::data->M.resize(run::mesh->elements->ndomains);
	run::data->f.resize(run::mesh->elements->ndomains);
	run::data->R.resize(run::mesh->elements->ndomains);
	run::data->primalSolution.resize(run::mesh->elements->ndomains);
	for (eslocal d = 0; d < run::mesh->elements->ndomains; d++) {
		run::data->K[d].rows = _domainDOFsSize[d];
		run::data->K[d].cols = _domainDOFsSize[d];
		run::data->f[d].resize(_domainDOFsSize[d]);
		run::data->R[d].resize(_domainDOFsSize[d]);
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = run::mesh->elements->domainDistribution[t]; d != run::mesh->elements->domainDistribution[t + 1]; ++d) {
			MatrixType mtype = _provider.getMatrixType(d);

			auto ebegin = run::mesh->elements->procNodes->cbegin() + run::mesh->elements->elementsDistribution[d];
			auto eend = run::mesh->elements->procNodes->cbegin() + run::mesh->elements->elementsDistribution[d + 1];

			eslocal Ksize = 0, RHSsize = 0;
			for (auto e = ebegin; e != eend; ++e) {
				RHSsize += e->size() * _DOFs;
				Ksize += getMatrixSize(e->size() * _DOFs, mtype);
			}

			for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
				if (
						run::mesh->boundaryRegions[r]->dimension &&
						run::mesh->boundaryRegions[r]->eintervalsDistribution[d] < run::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]) {

					eslocal begin = run::mesh->boundaryRegions[r]->eintervals[run::mesh->boundaryRegions[r]->eintervalsDistribution[d]].begin;
					eslocal end = run::mesh->boundaryRegions[r]->eintervals[run::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					auto enodes = run::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

					for (eslocal i = begin; i < end; ++i, ++enodes) {
						RHSsize += enodes->size() * _DOFs;
						Ksize += getMatrixSize(enodes->size() * _DOFs, mtype);;
					}
				}
			}

			std::vector<eslocal> permK(Ksize), permRHS(RHSsize);
			std::vector<IJ> KPattern(Ksize);
			std::vector<eslocal> RHSPattern(RHSsize), ROW, COL;

			IJ *Koffset = KPattern.data();
			eslocal *RHSoffset = RHSPattern.data();

			auto insert = [&] (serializededata<eslocal, eslocal>::const_iterator &enodes) {
				eslocal *_RHS = RHSoffset;
				for (auto n = enodes->begin(); n != enodes->end(); ++n, ++RHSoffset) {
					auto DOFs = (_DOFMap->begin() + *n)->begin();
					while (*DOFs != d + run::mesh->elements->firstDomain) {
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

			for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
				if (
						run::mesh->boundaryRegions[r]->dimension &&
						run::mesh->boundaryRegions[r]->eintervalsDistribution[d] < run::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]) {

					eslocal begin = run::mesh->boundaryRegions[r]->eintervals[run::mesh->boundaryRegions[r]->eintervalsDistribution[d]].begin;
					eslocal end = run::mesh->boundaryRegions[r]->eintervals[run::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					auto enodes = run::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

					for (eslocal i = begin; i < end; ++i, ++enodes) {
						insert(enodes);
						Koffset += getMatrixSize(enodes->size() * _DOFs, mtype);
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

			ROW.reserve(_domainDOFsSize[d] + 1);
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

			_KPermutation[d].swap(permK);
			_RHSPermutation[d].swap(permRHS);

			ROW.push_back(COL.size() + 1);

			run::data->K[d].nnz = COL.size();
			run::data->K[d].mtype = _provider.getMatrixType(d);
			switch (run::data->K[d].mtype) {
			case MatrixType::REAL_UNSYMMETRIC: run::data->K[d].type = 'G'; break;
			default: run::data->K[d].type = 'S';
			}
			run::data->K[d].CSR_V_values.resize(COL.size());
			run::data->K[d].CSR_I_row_indices.swap(ROW);
			run::data->K[d].CSR_J_col_indices.swap(COL);
			run::data->M[d] = run::data->K[d];
			run::data->M[d].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
			run::data->M[d].type = 'S';
		}
	}

//	std::cout << _KPermutation.front();
//	std::cout << _RHSPermutation.front();
}

void UniformNodesFETIComposer::buildB1Pattern()
{
	eslocal doffset = 0;

	run::data->B1.resize(run::mesh->elements->ndomains);
	run::data->B1c.resize(run::mesh->elements->ndomains);
	run::data->B1subdomainsMap.resize(run::mesh->elements->ndomains);
	run::data->B1duplicity.resize(run::mesh->elements->ndomains);
	run::data->LB.resize(run::mesh->elements->ndomains);
	run::data->block.resize(3);

	auto dmap = _DOFMap->begin();
	for (size_t i = 0, prev = 0; i < _dirichletMap.size(); prev = _dirichletMap[i++] / _DOFs) {
		dmap += (_dirichletMap[i] / _DOFs) - prev;
		for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
			if (run::mesh->elements->firstDomain <= *d && *d < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
				++doffset;
			}
			if (!_configuration.redundant_lagrange) {
				break;
			}
		}
	}

	eslocal dsize = Communication::exscan(doffset);

	dmap = _DOFMap->begin();
	for (size_t i = 0, prev = 0; i < _dirichletMap.size(); prev = _dirichletMap[i++] / _DOFs) {
		dmap += (_dirichletMap[i] / _DOFs) - prev;
		for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
			if (run::mesh->elements->firstDomain <= *d && *d < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
				run::data->B1[*d - run::mesh->elements->firstDomain].I_row_indices.push_back(doffset + 1);
				run::data->B1[*d - run::mesh->elements->firstDomain].J_col_indices.push_back(*(d + _dirichletMap[i] % _DOFs + 1) + 1);
				run::data->B1clustersMap.push_back({ doffset, environment->MPIrank });
				++doffset;
			}
			if (!_configuration.redundant_lagrange) {
				break;
			}
		}
	}

	_domainDirichletSize.resize(run::mesh->elements->ndomains);
	for (eslocal d = 0; d < run::mesh->elements->ndomains; ++d) {
		run::data->B1[d].V_values.resize(run::data->B1[d].I_row_indices.size(), 1);
		run::data->B1c[d].resize(run::data->B1[d].I_row_indices.size());
		run::data->B1duplicity[d].resize(run::data->B1[d].I_row_indices.size(), 1);
		_domainDirichletSize[d] = run::data->B1[d].I_row_indices.size();
	}

	run::data->block[DataHolder::CONSTRAINT::DIRICHLET] = dsize;

	eslocal goffset = 0;

	dmap = _DOFMap->begin();
	auto nranks = run::mesh->nodes->ranks->begin();
	auto exclude = _dirichletMap.begin();
	std::vector<std::vector<eslocal> > sBuffer(run::mesh->neighbours.size()), rBuffer(run::mesh->neighbours.size());

	auto send = [&] () {
		eslocal noffset = 0;
		for (auto r = nranks->begin() + 1; r != nranks->end(); ++r) {
			while (run::mesh->neighbours[noffset] < *r) {
				++noffset;
			}
			sBuffer[noffset].push_back(goffset);
		}
	};

	for (size_t n = 0; n < run::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (dmap->size() / (1 + _DOFs) > 1) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				eslocal ndomains = dmap->size() / (1 + _DOFs);
				if (_configuration.redundant_lagrange) {
					while (exclude != _dirichletMap.end() && *exclude < n * _DOFs + dof) {
						++exclude;
					}
					if (exclude == _dirichletMap.end() || n * _DOFs + dof != *exclude) {
						if (*nranks->begin() == environment->MPIrank) {
							send();
							goffset += ndomains * (ndomains - 1) / 2;
						}
					}
				} else {
					if (*nranks->begin() == environment->MPIrank) {
						send();
						goffset += ndomains - 1;
					}
				}
			}
		}
	}

	eslocal gsize = Communication::exscan(goffset);

	for (size_t n = 0; n < run::mesh->neighbours.size(); ++n) {
		for (size_t i = 0; i < sBuffer[n].size(); ++i) {
			sBuffer[n][i] += goffset;
		}
	}

	if (!Communication::receiveLowerUnknownSize(sBuffer, rBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange gluing offsets";
	}

	std::vector<eslocal> dDistribution = run::mesh->elements->gatherDomainsProcDistribution();

	dmap = _DOFMap->begin();
	nranks = run::mesh->nodes->ranks->begin();
	exclude = _dirichletMap.begin();

	auto push = [&] (eslocal d, eslocal lambda, eslocal dof, double value, eslocal domains) {
		d -= run::mesh->elements->firstDomain;
		run::data->B1[d].I_row_indices.push_back(lambda + 1);
		run::data->B1[d].J_col_indices.push_back(dof + 1);
		run::data->B1[d].V_values.push_back(value);
		run::data->B1duplicity[d].push_back(1. / domains);
	};

	auto fill = [&] (eslocal d1, eslocal d2, eslocal lambda, eslocal dof1, eslocal dof2, eslocal domains) {
		if (run::mesh->elements->firstDomain <= d1 && d1 < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
			push(d1, lambda, dof1, 1, domains);
		}
		if (run::mesh->elements->firstDomain <= d2 && d2 < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
			push(d2, lambda, dof2, -1, domains);
		}
		if (
				(run::mesh->elements->firstDomain <= d1 && d1 < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) ||
				(run::mesh->elements->firstDomain <= d2 && d2 < run::mesh->elements->firstDomain + run::mesh->elements->ndomains)) {

			run::data->B1clustersMap.push_back({ lambda, environment->MPIrank });

			if (d1 <run::mesh->elements->firstDomain || run::mesh->elements->firstDomain + run::mesh->elements->ndomains <= d1) {
				run::data->B1clustersMap.back().push_back(std::lower_bound(dDistribution.begin(), dDistribution.end(), d1 + 1) - dDistribution.begin() - 1);
			}
			if (d2 <run::mesh->elements->firstDomain || run::mesh->elements->firstDomain + run::mesh->elements->ndomains <= d2) {
				run::data->B1clustersMap.back().push_back(std::lower_bound(dDistribution.begin(), dDistribution.end(), d2 + 1) - dDistribution.begin() - 1);
			}
		}
	};

	std::vector<eslocal> roffset(run::mesh->neighbours.size());
	for (size_t n = 0; n < run::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (dmap->size() / (1 + _DOFs) > 1) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				if (_configuration.redundant_lagrange) {
					while (exclude != _dirichletMap.end() && *exclude < n * _DOFs + dof) {
						++exclude;
					}
				}
				eslocal lambda = goffset;
				eslocal ndomains = dmap->size() / (1 + _DOFs);
				if (*nranks->begin() != environment->MPIrank) {
					if (!_configuration.redundant_lagrange || exclude == _dirichletMap.end() || n * _DOFs + dof != *exclude) {
						eslocal noffset = 0;
						while (run::mesh->neighbours[noffset] < *nranks->begin()) {
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
						if (*nranks->begin() == environment->MPIrank) {
							goffset += ndomains * (ndomains - 1) / 2;
						}
					}
				} else {
					for (auto d1 = dmap->begin(), d2 = dmap->begin() + 1 + _DOFs; d2 != dmap->end(); d1 = d2, d2 += 1 + _DOFs, ++lambda) {
						fill(*d1, *d2, dsize + lambda, *(d1 + 1 + dof), *(d2 + 1 + dof), ndomains);
					}
					if (*nranks->begin() == environment->MPIrank) {
						goffset += ndomains - 1;
					}
				}
			}
		}
	}


	for (eslocal d = 0; d < run::mesh->elements->ndomains; ++d) {
		run::data->B1c[d].resize(run::data->B1[d].I_row_indices.size());

		run::data->B1[d].nnz = run::data->B1[d].I_row_indices.size();
		run::data->B1[d].rows = dsize + gsize;
		run::data->B1[d].cols = _domainDOFsSize[d];
		run::data->LB[d].resize(run::data->B1[d].nnz, -std::numeric_limits<double>::infinity());

		run::data->B1subdomainsMap[d] = run::data->B1[d].I_row_indices;
		for (size_t i = 0; i < run::data->B1subdomainsMap[d].size(); ++i) {
			--run::data->B1subdomainsMap[d][i];
		}
	}

	run::data->block[DataHolder::CONSTRAINT::EQUALITY_CONSTRAINTS] = dsize + gsize;
}

void UniformNodesFETIComposer::buildB0Pattern()
{
	run::data->B0.resize(run::mesh->elements->ndomains);

	for (eslocal d = 0; d < run::mesh->elements->ndomains; ++d) {
		run::data->B0[d].rows = 0;
		run::data->B0[d].cols = 0;
		run::data->B0[d].nnz = 0;
	}
}

void UniformNodesFETIComposer::assemble(Matrices matrices, const SolverParameters &parameters)
{
	_controler.nextTime();

	#pragma omp parallel for
	for  (eslocal d = 0; d < run::mesh->elements->ndomains; d++) {

		size_t KIndex = 0, RHSIndex = 0;
		double KReduction = 1, RHSReduction = 1; //_step.internalForceReduction;
		Controler::InstanceFiller filler;

		switch (_provider.getMatrixType(d)) {
		case MatrixType::REAL_UNSYMMETRIC:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						run::data->f[d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						run::data->R[d][_RHSPermutation[d][RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = 0; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							run::data->K[d].CSR_V_values[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							run::data->M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						run::data->f[d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						run::data->R[d][_RHSPermutation[d][RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							run::data->K[d].CSR_V_values[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							run::data->M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		}

		clearMatrices(matrices, d);
		filler.begin = run::mesh->elements->elementsDistribution[d];
		filler.end = run::mesh->elements->elementsDistribution[d + 1];

		_controler.processElements(matrices, parameters, filler);

		KReduction = 1; //_step.internalForceReduction;

		for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
			if (run::mesh->boundaryRegions[r]->distribution.size()) {
				if (run::mesh->boundaryRegions[r]->eintervalsDistribution[d] < run::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]) {
					filler.begin = run::mesh->boundaryRegions[r]->eintervals[run::mesh->boundaryRegions[r]->eintervalsDistribution[d]].begin;
					filler.end = run::mesh->boundaryRegions[r]->eintervals[run::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					_controler.processBoundary(matrices, parameters, r, filler);
				}
			}
		}
	}
}

void UniformNodesFETIComposer::setDirichlet()
{
	std::vector<double> values(_dirichletMap.size());
	_controler.dirichletValues(values);

	std::vector<eslocal> doffset(run::mesh->elements->ndomains);

	auto dmap = _DOFMap->begin();
	for (size_t i = 0, prev = 0; i < _dirichletMap.size(); prev = _dirichletMap[i++] / _DOFs) {
		dmap += (_dirichletMap[i] / _DOFs) - prev;
		for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
			if (run::mesh->elements->firstDomain <= *d && *d < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
				run::data->B1c[*d - run::mesh->elements->firstDomain][doffset[*d - run::mesh->elements->firstDomain]++] = values[_dirichletPermutation[i]];
			}
			if (!_configuration.redundant_lagrange) {
				break;
			}
		}
	}
}

void UniformNodesFETIComposer::synchronize()
{
	if (_configuration.scaling && _configuration.redundant_lagrange) {
		updateDuplicity();
	}
}

void UniformNodesFETIComposer::updateDuplicity()
{
	std::vector<std::vector<double> > diagonals(run::mesh->elements->ndomains);
	#pragma omp parallel for
	for  (eslocal d = 0; d < run::mesh->elements->ndomains; d++) {
		diagonals[d] = run::data->K[d].getDiagonal();
	}

	std::vector<std::vector<double> > sBuffer(run::mesh->neighbours.size()), rBuffer(run::mesh->neighbours.size());

	auto dmap = _DOFMap->begin();
	auto nranks = run::mesh->nodes->ranks->begin();
	std::vector<double> buffer;

	for (eslocal n = 0; n < run::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (nranks->size() > 1) {
			buffer.clear();
			for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
				if (run::mesh->elements->firstDomain <= *d && *d < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
					for (int dof = 0; dof < _DOFs; ++dof) {
						buffer.push_back(diagonals[*d - run::mesh->elements->firstDomain][*(d + 1 + dof)]);
					}
				}
			}

			eslocal noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != environment->MPIrank) {
					while (run::mesh->neighbours[noffset] < *r) {
						++noffset;
					}
					sBuffer[noffset].insert(sBuffer[noffset].end(), buffer.begin(), buffer.end());
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange diagonal values";
	}

	dmap = _DOFMap->begin();
	nranks = run::mesh->nodes->ranks->begin();
	std::vector<eslocal> roffset(run::mesh->neighbours.size());
	std::vector<eslocal> dDistribution = run::mesh->elements->gatherDomainsProcDistribution();
	auto exclude = _dirichletMap.begin();

	auto push = [&] (eslocal d1, eslocal d2, double v1, double v2, double sum) {
		if (run::mesh->elements->firstDomain <= d1 && d1 < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
			run::data->B1duplicity[d1 - run::mesh->elements->firstDomain].push_back(v2 / sum);
		}
		if (run::mesh->elements->firstDomain <= d2 && d2 < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
			run::data->B1duplicity[d2 - run::mesh->elements->firstDomain].push_back(v1 / sum);
		}
	};

	for (size_t d = 0; d < run::data->B1duplicity.size(); d++) {
		run::data->B1duplicity[d].resize(_domainDirichletSize[d]);
	}
	for (eslocal n = 0; n < run::mesh->nodes->size; ++n, ++dmap, ++nranks) {
		if (dmap->size() / (1 + _DOFs) > 1) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				buffer.clear();

				for (auto d = dmap->begin(); d != dmap->end(); d += 1 + _DOFs) {
					if (run::mesh->elements->firstDomain <= *d && *d < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
						buffer.push_back(diagonals[*d - run::mesh->elements->firstDomain][*(d + 1 + dof)]);
					} else {
						eslocal r = std::lower_bound(dDistribution.begin(), dDistribution.end(), *d + 1) - dDistribution.begin() - 1;
						eslocal noffset = 0;
						while (run::mesh->neighbours[noffset] < r) {
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
					eslocal v1 = 0, v2 = 1;
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
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<double> &solution = _controler.getSolutionStore();

	std::vector<std::vector<std::vector<double> > > sBuffer(threads);
	std::vector<std::vector<double> > rBuffer(run::mesh->neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t i = run::mesh->nodes->distribution[t];
		auto nranks = run::mesh->nodes->ranks->begin(t);

		std::vector<std::vector<double> > tBuffer(run::mesh->neighbours.size());

		for (auto map = _DOFMap->begin(t); map != _DOFMap->end(t); ++map, ++i, ++nranks) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				solution[i * _DOFs + dof] = 0;
			}
			for (auto d = map->begin(); d != map->end(); d += 1 + _DOFs) {
				if (run::mesh->elements->firstDomain <= *d && *d < run::mesh->elements->firstDomain + run::mesh->elements->ndomains) {
					for (int dof = 0; dof < _DOFs; ++dof) {
						solution[i * _DOFs + dof] += run::data->primalSolution[*d - run::mesh->elements->firstDomain][*(d + 1 + dof)] / ((map->end() - map->begin()) / (1 + _DOFs));
					}
				}
			}

			eslocal noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != environment->MPIrank) {
					while (run::mesh->neighbours[noffset] < *r) {
						++noffset;
					}

					for (size_t dof = 0; dof < _DOFs; ++dof) {
						tBuffer[noffset].push_back(solution[i * _DOFs + dof]);
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

	if (!Communication::exchangeKnownSize(sBuffer[0], rBuffer, run::mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: synchronize solution.";
	}

	std::vector<eslocal> roffset(run::mesh->neighbours.size());
	auto nranks = run::mesh->nodes->ranks->begin();
	for (eslocal n = 0; n < run::mesh->nodes->size; ++n, ++nranks) {
		eslocal noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			if (*r != environment->MPIrank) {
				while (run::mesh->neighbours[noffset] < *r) {
					++noffset;
				}

				for (size_t dof = 0; dof < _DOFs; ++dof) {
					solution[n * _DOFs + dof] += rBuffer[noffset][roffset[noffset]];
				}
				++roffset[noffset];
			}
		}
	}
}
