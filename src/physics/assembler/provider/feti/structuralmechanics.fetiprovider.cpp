
#include "esinfo/meshinfo.h"
#include "physics/assembler/dataholder.h"
#include "structuralmechanics.fetiprovider.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/matrixtype.h"
#include "config/ecf/physics/structuralmechanics.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/fetidatastore.h"

#include "solver/generic/SparseMatrix.h"
#include "solver/specific/sparsesolvers.h"

#include <algorithm>

using namespace espreso;

StructuralMechanicsFETIProvider::StructuralMechanicsFETIProvider(DataHolder *data, StructuralMechanicsLoadStepConfiguration &configuration)
: FETIProvider(data, configuration), _configuration(configuration)
{
	if (configuration.feti.regularization == FETI_REGULARIZATION::ANALYTIC) {
		prepareRegularization();
	}
}

MatrixType StructuralMechanicsFETIProvider::getMatrixType() const
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

MatrixType StructuralMechanicsFETIProvider::getMatrixType(esint domain) const
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void StructuralMechanicsFETIProvider::prepareRegularization()
{
	size_t clusters = *std::max_element(info::mesh->elements->clusters.begin(), info::mesh->elements->clusters.end()) + 1;

	_cCenter = _cNorm = std::vector<Point>(clusters, Point(0, 0, 0));
	_cr44 = _cr45 = _cr46 = _cr55 = _cr56 = std::vector<double>(clusters, 0);
	_cNp = std::vector<size_t>(clusters, 0);

	_dCenter = _dNorm = std::vector<Point>(info::mesh->elements->ndomains, Point(0, 0, 0));
	_dr44 = _dr45 = _dr46 = _dr55 = _dr56 = std::vector<double>(info::mesh->elements->ndomains, 0);
	_dNp = std::vector<size_t>(info::mesh->elements->ndomains, 0);

	std::vector<double> cbuffer1(info::mesh->elements->ndomains, 0), cbuffer2(info::mesh->elements->ndomains, 0), cbuffer3(info::mesh->elements->ndomains, 0);

	// Get center
	#pragma omp parallel for
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		Point center;
		for (size_t i = 0; i < info::mesh->nodes->dintervals[p].size(); i++) {
			for (esint n = info::mesh->nodes->dintervals[p][i].begin; n < info::mesh->nodes->dintervals[p][i].end; ++n) {
				center += info::mesh->nodes->coordinates->datatarray()[n];
			}
		}
		_dCenter[p] = center;
	}

	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cCenter[info::mesh->elements->clusters[p]] += _dCenter[p];
		_dNp[p] = info::mesh->nodes->dintervals[p].back().DOFOffset + info::mesh->nodes->dintervals[p].back().end - info::mesh->nodes->dintervals[p].back().begin;
		_dCenter[p] = _dCenter[p] / _dNp[p];
		_cNp[info::mesh->elements->clusters[p]] += _dNp[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cCenter[c] /= _cNp[c];
	}

	// Compute norm of column 4 (norm.x)
	#pragma omp parallel for
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		double pnorm = 0, pcnorm = 0;
		for (size_t i = 0; i < info::mesh->nodes->dintervals[p].size(); i++) {
			for (esint n = info::mesh->nodes->dintervals[p][i].begin; n < info::mesh->nodes->dintervals[p][i].end; ++n) {
				Point dp = info::mesh->nodes->coordinates->datatarray()[n] - _dCenter[p];
				pnorm += dp.x * dp.x + dp.y * dp.y;
				Point cp = info::mesh->nodes->coordinates->datatarray()[n] - _cCenter[info::mesh->elements->clusters[p]];
				pcnorm += cp.x * cp.x + cp.y * cp.y;
			}
		}
		_dNorm[p].x = std::sqrt(pnorm);
		cbuffer1[p] += pcnorm;
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cNorm[info::mesh->elements->clusters[p]].x += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].x = std::sqrt(_cNorm[c].x);
	}

	// Compute coefficient r44, r45
	cbuffer1 = cbuffer2 = std::vector<double>(info::mesh->elements->ndomains, 0);
	#pragma omp parallel for
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		size_t c = info::mesh->elements->clusters[p];
		for (size_t i = 0; i < info::mesh->nodes->dintervals[p].size(); i++) {
			for (esint n = info::mesh->nodes->dintervals[p][i].begin; n < info::mesh->nodes->dintervals[p][i].end; ++n) {
				Point dp = info::mesh->nodes->coordinates->datatarray()[n] - _dCenter[p];
				_dr44[p] += (-dp.y / _dNorm[p].x) * (-dp.y / _dNorm[p].x) + (dp.x / _dNorm[p].x) * (dp.x / _dNorm[p].x);
				_dr45[p] += (-dp.y / _dNorm[p].x) * (-dp.z);

				Point cp = info::mesh->nodes->coordinates->datatarray()[n] - _cCenter[c];
				cbuffer1[p] += (-cp.y / _cNorm[c].x) * (-cp.y / _cNorm[c].x) + (cp.x / _cNorm[c].x) * (cp.x / _cNorm[c].x);
				cbuffer2[p] += (-cp.y / _cNorm[c].x) * (-cp.z);
			}
		}
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cr44[info::mesh->elements->clusters[p]] += cbuffer1[p];
		_cr45[info::mesh->elements->clusters[p]] += cbuffer2[p];
	}

	// Compute norm of column 5 (norm.y)
	cbuffer1 = std::vector<double>(info::mesh->elements->ndomains, 0);
	#pragma omp parallel for
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		double dnorm = 0, cnorm = 0;
		size_t c = info::mesh->elements->clusters[p];
		for (size_t i = 0; i < info::mesh->nodes->dintervals[p].size(); i++) {
			for (esint n = info::mesh->nodes->dintervals[p][i].begin; n < info::mesh->nodes->dintervals[p][i].end; ++n) {
				Point dp = info::mesh->nodes->coordinates->datatarray()[n] - _dCenter[p];
				dnorm += (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x)) * (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x));
				dnorm += (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x)) * (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x));
				dnorm += dp.x * dp.x;

				Point cp = info::mesh->nodes->coordinates->datatarray()[n] - _cCenter[c];
				cnorm += (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x)) * (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x));
				cnorm += (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x)) * (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x));
				cnorm += cp.x * cp.x;
			}
		}
		_dNorm[p].y = std::sqrt(dnorm);
		cbuffer1[p] = cnorm;
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cNorm[info::mesh->elements->clusters[p]].y += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].y = std::sqrt(_cNorm[c].y);
	}

	// Compute coefficient r46, r55, r56
	cbuffer1 = cbuffer2 = cbuffer3 = std::vector<double>(info::mesh->elements->ndomains, 0);
	#pragma omp parallel for
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		double c5;
		size_t c = info::mesh->elements->clusters[p];
		for (size_t i = 0; i < info::mesh->nodes->dintervals[p].size(); i++) {
			for (esint n = info::mesh->nodes->dintervals[p][i].begin; n < info::mesh->nodes->dintervals[p][i].end; ++n) {
				Point dp = info::mesh->nodes->coordinates->datatarray()[n] - _dCenter[p];
				_dr46[p] += (dp.x / _dNorm[p].x) * (-dp.z);
				c5 = (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x)) / _dNorm[p].y;
				_dr55[p] += c5 * c5;
				_dr56[p] += c5 * 0;
				c5 = (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x)) / _dNorm[p].y;
				_dr55[p] += c5 * c5;
				_dr56[p] += c5 * (-dp.z);
				c5 = ( dp.x -                                           0) / _dNorm[p].y;
				_dr55[p] += c5 * c5;
				_dr56[p] += c5 * dp.y;

				Point cp = info::mesh->nodes->coordinates->datatarray()[n] - _cCenter[c];
				cbuffer1[p] += (cp.x / _cNorm[c].x) * (-cp.z);
				c5 = (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x)) / _cNorm[c].y;
				cbuffer2[p] += c5 * c5;
				cbuffer3[p] += c5 * 0;
				c5 = (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x)) / _cNorm[c].y;
				cbuffer2[p] += c5 * c5;
				cbuffer3[p] += c5 * (-cp.z);
				c5 = ( cp.x -                                           0) / _cNorm[c].y;
				cbuffer2[p] += c5 * c5;
				cbuffer3[p] += c5 * cp.y;
			}
		}
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cr46[info::mesh->elements->clusters[p]] += cbuffer1[p];
		_cr55[info::mesh->elements->clusters[p]] += cbuffer2[p];
		_cr56[info::mesh->elements->clusters[p]] += cbuffer3[p];
	}

	// Compute norm of column 6 (norm.z)
	cbuffer1 = std::vector<double>(info::mesh->elements->ndomains, 0);
	#pragma omp parallel for
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		double dnorm = 0, cnorm = 0, c6;
		size_t c = info::mesh->elements->clusters[p];
		for (size_t i = 0; i < info::mesh->nodes->dintervals[p].size(); i++) {
			for (esint n = info::mesh->nodes->dintervals[p][i].begin; n < info::mesh->nodes->dintervals[p][i].end; ++n) {
				Point dp = info::mesh->nodes->coordinates->datatarray()[n] - _dCenter[p];
				c6 =     0 - _dr56[p] / _dr55[p] * (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x)) / _dNorm[p].y - _dr46[p] / _dr44[p] * (-dp.y / _dNorm[p].x);
				dnorm += c6 * c6;
				c6 = -dp.z - _dr56[p] / _dr55[p] * (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x)) / _dNorm[p].y - _dr46[p] / _dr44[p] * ( dp.x / _dNorm[p].x);
				dnorm += c6 * c6;
				c6 =  dp.y - _dr56[p] / _dr55[p] * ( dp.x -                                           0) / _dNorm[p].y - _dr46[p] / _dr44[p] * (    0 / _dNorm[p].x);
				dnorm += c6 * c6;

				Point cp = info::mesh->nodes->coordinates->datatarray()[n] - _cCenter[c];
				c6 =     0 - _cr56[c] / _cr55[c] * (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x)) / _cNorm[c].y - _cr46[c] / _cr44[c] * (-cp.y / _cNorm[c].x);
				cnorm += c6 * c6;
				c6 = -cp.z - _cr56[c] / _cr55[c] * (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x)) / _cNorm[c].y - _cr46[c] / _cr44[c] * ( cp.x / _cNorm[c].x);
				cnorm += c6 * c6;
				c6 =  cp.y - _cr56[c] / _cr55[c] * ( cp.x -                                           0) / _cNorm[c].y - _cr46[c] / _cr44[c] * (    0 / _cNorm[c].x);
				cnorm += c6 * c6;
			}
		}
		_dNorm[p].z = std::sqrt(dnorm);
		cbuffer1[p] = cnorm;
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cNorm[info::mesh->elements->clusters[p]].z += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].z = std::sqrt(_cNorm[c].z);
	}
}





