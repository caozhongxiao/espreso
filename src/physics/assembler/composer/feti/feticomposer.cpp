
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "config/ecf/physics/physics.h"
#include "basis/containers/serializededata.h"
#include "physics/assembler/dataholder.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/assembler/assembler.h"
#include "physics/assembler/provider/feti/fetiprovider.h"
#include "feticomposer.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/nodestore.h"
#include "solver/generic/SparseMatrix.h"

#include "wrappers/bem/bemwrapper.h"

using namespace espreso;

FETIComposer::FETIComposer(Controller &controler, FETIProvider &provider, FETISolverConfiguration &configuration)
: Composer(controler), _provider(provider), _configuration(configuration)
{
	_BEMDomain.resize(info::mesh->elements->ndomains);
	if (BEM4I::isLinked()) {
		for (esint d = 0; d < info::mesh->elements->ndomains; ++d) {
			_BEMDomain[d] = isBEMDomain(d);
		}
	}
}

void FETIComposer::assemble(Matrices matrices, const SolverParameters &parameters)
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
						if ((matrices & Matrices::M) && filler.Me.rows() && r / filler.Me.rows() == c / filler.Me.rows()) {
							data->M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r % filler.Me.rows(), c % filler.Me.rows());;
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
						if ((matrices & Matrices::M) && filler.Me.rows() && r / filler.Me.rows() == c / filler.Me.rows()) {
							data->M[d].CSR_V_values[_KPermutation[d][KIndex]] += filler.Me(r % filler.Me.rows(), c % filler.Me.rows());;
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

void FETIComposer::fillSolution()
{
	for (size_t i = 0; i < _BEMDomain.size(); i++) {
		if (_BEMDomain[i]) {
			data->primalSolution[i].resize(_domainDOFsSize[i]);
			_controler.fillBEMInterior(i, data->primalSolution[i].data());
		}
	}
	avgGather(_controler.solution()->data, data->primalSolution);
}

void FETIComposer::KplusAlfaM(double alfa)
{
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		data->K[d].MatAddInPlace(data->M[d], 'N', alfa);
	}
}

void FETIComposer::alfaKplusBetaM(double alfa, double beta)
{
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		data->K[d].MatScale(alfa);
		data->K[d].MatAddInPlace(data->M[d], 'N', alfa);
	}
}

bool FETIComposer::isBEMDomain(esint domain)
{
	auto eregions = (info::mesh->elements->regions->begin() + info::mesh->elements->elementsDistribution[domain])->begin();
	for (int byte = 0; byte < info::mesh->elements->regionMaskSize; ++byte) {
		for (size_t bit = 0; bit < sizeof(esint); bit++) {
			if (eregions[byte] & 1 << bit) {
				auto region = info::mesh->elementsRegions[byte * sizeof(esint) + bit];
				if (_controler.configuration().discretization.find(region->name) != _controler.configuration().discretization.end()) {
					return true;
				}
			}
		}
	}
	return false;
}


void FETIComposer::apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x)
{
	std::vector<std::vector<double> > _res(matrices.size()), _x(matrices.size());

	for (size_t i = 0; i < matrices.size(); i++) {
		_res[i].resize(matrices[i].rows);
		_x[i].resize(matrices[i].rows);
	}

	duply(x, _x);

	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		matrices[d].MatVec(_x[d], _res[d], 'N', 0, 0, 0);
	}

	gather(result, _res);
}

void FETIComposer::enrichRHS(double alfa, NodeData* x)
{
	std::vector<std::vector<double> > _x(data->f.size());

	for (size_t i = 0; i < data->f.size(); i++) {
		_x[i].resize(data->f[i].size());
	}
	divide(x->data, _x);

	#pragma omp parallel for
	for (size_t d = 0; d < data->f.size(); d++) {
		for (size_t i = 0; i < data->f[d].size(); ++i) {
			data->f[d][i] += _x[d][i];
		}
	}
}

double FETIComposer::residualNormNumerator()
{
	double square = 0;
	std::vector<double> f, btlambda;
	gather(f, data->f);
	gather(btlambda, data->dualSolution);
	for (size_t i = _foreignDOFs; i < f.size(); i++) {
		square += (f[i] - btlambda[i]) * (f[i] - btlambda[i]);
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}

double FETIComposer::residualNormDenominator()
{
	double square = 0;
	std::vector<double> f, r;
	gather(f, data->origF);
	gather(r, data->R);
	for (size_t i = _foreignDOFs, d = 0; i < f.size(); i++) {
		while (d < _dirichletMap.size() && (size_t)_dirichletMap[d] < i) { d++; }
		if (d == _dirichletMap.size() || (size_t)_dirichletMap[d] != i) {
			square += f[i] * f[i];
		} else {
			square += r[i] * r[i];
		}
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}
