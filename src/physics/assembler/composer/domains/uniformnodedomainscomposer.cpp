
#include "uniformnodedomainscomposer.h"

#include "../../controlers/controler.h"

#include "../../../instance.h"
#include "../../../step.h"

#include "../../../../basis/containers/serializededata.h"
#include "../../../../basis/matrices/denseMatrix.h"
#include "../../../../basis/logging/logging.h"
#include "../../../../basis/utilities/utils.h"
#include "../../../../basis/utilities/communication.h"

#include "../../../../config/ecf/environment.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"
#include "../../../../mesh/store/nodestore.h"

#include "../../../../solver/generic/SparseMatrix.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void UniformNodeDomainsComposer::initDOFs()
{
	size_t threads = environment->OMP_NUM_THREADS;

	// nID, domain
	std::vector<std::vector<std::pair<eslocal, eslocal> > > ntodomains(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<eslocal, eslocal> > tdata;

		for (eslocal d = _mesh.elements->domainDistribution[t]; d != _mesh.elements->domainDistribution[t + 1]; ++d) {
			std::vector<eslocal> dnodes(
					(_mesh.elements->procNodes->begin() + _mesh.elements->elementsDistribution[d])->begin(),
					(_mesh.elements->procNodes->begin() + _mesh.elements->elementsDistribution[d + 1])->begin());

			Esutils::sortAndRemoveDuplicity(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<eslocal, eslocal>(dnodes[i], d));
			}
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
		std::vector<eslocal> tDOFs(_mesh.elements->ndomains);

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			tDOFs[ntodomains[0][n].second] += _DOFs;
		}

		DOFs[t].swap(tDOFs);
	}

	Esutils::sizesToOffsets(DOFs);

	std::vector<std::vector<std::vector<eslocal> > > sBuffer(threads);
	std::vector<std::vector<eslocal> > rBuffer(_mesh.neighboursWithMe.size());


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nranks = _mesh.nodes->ranks->begin() + ntodomains[0][ndistribution[t]].first;

		std::vector<std::vector<eslocal> > tBuffer(_mesh.neighboursWithMe.size());

		size_t n = ndistribution[t];
		while (n < ndistribution[t + 1]) {
			size_t begin = n++;
			while (n < ntodomains[0].size() && ntodomains[0][n].first == ntodomains[0][n - 1].first) {
				++n;
			}

			eslocal noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (_mesh.neighboursWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(n - begin);
				for (size_t i = begin; i < n; i++) {
					tBuffer[noffset].push_back(_mesh.elements->firstDomain + ntodomains[0][i].second);
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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh.neighboursWithMe)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange uniform DOFs.";
	}

	// domain0, DOF0, DOF1, ..., DOFn, domain1, DOF0, DOF1, ..., DOFn, ...; domain0, ...
	std::vector<eslocal> DOFDistribution(1), DOFData;

	// TODO: make it parallel
	// parallelization is possible if node order will be kept as: boundary first!
	// now we prefer generality
	auto nranks = _mesh.nodes->ranks->begin();
	std::vector<eslocal> roffset(rBuffer.size());
	for (eslocal n = 0; n < _mesh.nodes->size; ++n, ++nranks) {
		eslocal noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (_mesh.neighboursWithMe[noffset] < *r) {
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

	std::vector<size_t> distribution = _mesh.nodes->distribution, datadistribution(threads + 1);
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

void UniformNodeDomainsComposer::initDirichlet()
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
	std::sort(_dirichletMap.begin(), _dirichletMap.end());
}

void UniformNodeDomainsComposer::buildPatterns()
{
	size_t threads = environment->OMP_NUM_THREADS;

	_KPermutation.resize(_mesh.elements->ndomains);
	_RHSPermutation.resize(_mesh.elements->ndomains);

	_instance.K.resize(_mesh.elements->ndomains);
	_instance.M.resize(_mesh.elements->ndomains);
	_instance.f.resize(_mesh.elements->ndomains);
	_instance.R.resize(_mesh.elements->ndomains);
	_instance.primalSolution.resize(_mesh.elements->ndomains);
	for (eslocal d = 0; d < _mesh.elements->ndomains; d++) {
		_instance.K[d].rows = _instance.domainDOFCount[d];
		_instance.K[d].cols = _instance.domainDOFCount[d];
		_instance.f[d].resize(_instance.domainDOFCount[d]);
		_instance.R[d].resize(_instance.domainDOFCount[d]);
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = _mesh.elements->domainDistribution[t]; d != _mesh.elements->domainDistribution[t + 1]; ++d) {
			MatrixType mtype = _controler.getMatrixType(d);

			auto ebegin = _mesh.elements->procNodes->cbegin() + _mesh.elements->elementsDistribution[d];
			auto eend = _mesh.elements->procNodes->cbegin() + _mesh.elements->elementsDistribution[d + 1];

			eslocal Ksize = 0, RHSsize = 0;
			for (auto e = ebegin; e != eend; ++e) {
				RHSsize += e->size() * _DOFs;
				Ksize += getMatrixSize(e->size() * _DOFs, mtype);
			}

			for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
				if (
						_mesh.boundaryRegions[r]->dimension &&
						_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {

					eslocal begin = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d]].begin;
					eslocal end = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					auto enodes = _mesh.boundaryRegions[r]->procNodes->cbegin() + begin;

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
					while (*DOFs != d + _mesh.elements->firstDomain) {
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

			for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
				if (
						_mesh.boundaryRegions[r]->dimension &&
						_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {

					eslocal begin = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d]].begin;
					eslocal end = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
					auto enodes = _mesh.boundaryRegions[r]->procNodes->cbegin() + begin;

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

			ROW.reserve(_instance.domainDOFCount[d] + 1);
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

			_instance.K[d].nnz = COL.size();
			_instance.K[d].mtype = _controler.getMatrixType(d);
			switch (_instance.K[d].mtype) {
			case MatrixType::REAL_UNSYMMETRIC: _instance.K[d].type = 'G'; break;
			default: _instance.K[d].type = 'S';
			}
			_instance.K[d].CSR_V_values.resize(COL.size());
			_instance.K[d].CSR_I_row_indices.swap(ROW);
			_instance.K[d].CSR_J_col_indices.swap(COL);
			_instance.M[d] = _instance.K[d];
			_instance.M[d].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
			_instance.M[d].type = 'S';
		}
	}

//	std::cout << _KPermutation.front();
//	std::cout << _RHSPermutation.front();
}

void UniformNodeDomainsComposer::assemble(Matrices matrices)
{
	_controler.nextTime();

	#pragma omp parallel for
	for  (size_t d = 0; d < _instance.domains; d++) {

		size_t KIndex = 0, RHSIndex = 0;
		double KReduction = 1, RHSReduction = _step.internalForceReduction;
		Controler::InstanceFiller filler;

		switch (_controler.getMatrixType(d)) {
		case MatrixType::REAL_UNSYMMETRIC:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						_instance.f[d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						_instance.R[d][_RHSPermutation[d][RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = 0; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							_instance.K[d].CSR_V_values[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							_instance.M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		default:
			filler.insert = [&] (size_t size) {
				for (size_t r = 0; r < size; ++r, ++RHSIndex) {
					if ((matrices & Matrices::f) && filler.fe.rows()) {
						_instance.f[d][_RHSPermutation[d][RHSIndex]] += RHSReduction * filler.fe(r, 0);
					}
					if ((matrices & Matrices::R) && filler.Re.rows()) {
						_instance.R[d][_RHSPermutation[d][RHSIndex]] += filler.Re(r, 0);
					}

					for (size_t c = r; c < size; ++c, ++KIndex) {
						if ((matrices & Matrices::K) && filler.Ke.rows()) {
							_instance.K[d].CSR_V_values[_KPermutation[d][KIndex]] += KReduction * filler.Ke(r, c);
						}
						if ((matrices & Matrices::M) && filler.Me.rows()) {
							_instance.M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r, c);
						}
					}
				}
			}; break;
		}

		clearMatrices(matrices, d);
		filler.begin = _mesh.elements->elementsDistribution[d];
		filler.end = _mesh.elements->elementsDistribution[d + 1];

		_controler.processElements(matrices, filler);

//		KReduction = _step.internalForceReduction;
//
//		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
//			if (_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {
//				filler.begin = _mesh.elements->elementsDistribution[d];
//				filler.end = _mesh.elements->elementsDistribution[d + 1];
//			}
//		}

//		auto boundary = [&] (Matrices restriction) {
//			for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
//				if (_mesh.boundaryRegions[r]->dimension == 2) {
//					if (_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {
//						eslocal begin = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d]].begin;
//						eslocal end = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
//						auto nodes = _mesh.boundaryRegions[r]->procNodes->cbegin() + begin;
//						for (eslocal i = begin; i < end; ++i, ++nodes) {
//							processFace(d, _mesh.boundaryRegions[r], restriction, i, Ke, Me, Re, fe);
//							switch (getMatrixType(domain)) {
//							case MatrixType::REAL_UNSYMMETRIC: fullInsert(nodes, _step.internalForceReduction); break;
//							default: upperInsert(nodes, _step.internalForceReduction);
//							}
//						}
//					}
//				}
//				if (_mesh.boundaryRegions[r]->dimension == 1) {
//					if (_mesh.boundaryRegions[r]->eintervalsDistribution[d] < _mesh.boundaryRegions[r]->eintervalsDistribution[d + 1]) {
//						eslocal begin = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d]].begin;
//						eslocal end = _mesh.boundaryRegions[r]->eintervals[_mesh.boundaryRegions[r]->eintervalsDistribution[d + 1] - 1].end;
//						auto nodes = _mesh.boundaryRegions[r]->procNodes->cbegin() + begin;
//						for (eslocal i = begin; i < end; ++i, ++nodes) {
//							processEdge(d, _mesh.boundaryRegions[r], restriction, i, Ke, Me, Re, fe);
//							switch (getMatrixType(domain)) {
//							case MatrixType::REAL_UNSYMMETRIC: fullInsert(nodes, _step.internalForceReduction); break;
//							default: upperInsert(nodes, _step.internalForceReduction);
//							}
//						}
//					}
//				}
//			}
//		};
	}
}

void UniformNodeDomainsComposer::setDirichlet()
{

}

void UniformNodeDomainsComposer::synchronize()
{

}

void UniformNodeDomainsComposer::fillSolution()
{
	std::vector<double> &solution = _controler.getSolutionStore();
}

