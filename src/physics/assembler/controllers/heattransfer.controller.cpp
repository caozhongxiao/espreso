
#include "heattransfer2d.controller.h"

#include "../../../globals/run.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../basis/matrices/matrixtype.h"
#include "../../../config/ecf/physics/heattransfer.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"
#include "../../../mesh/store/boundaryregionstore.h"

#include "../../../basis/utilities/communication.h"
#include "../../../basis/utilities/utils.h"
#include "../../../globals/time.h"

using namespace espreso;

MatrixType HeatTransferControler::getMatrixType(size_t domain) const
{
//	if (_step.tangentMatrixCorrection) {
//		return MatrixType::REAL_UNSYMMETRIC;
//	}

	if (_configuration.translation_motions.size()) {
		for (auto it = _configuration.translation_motions.begin(); it != _configuration.translation_motions.end(); ++it) {
			ElementsRegionStore *region = run::mesh->eregion(it->first);
			for (eslocal i = run::mesh->elements->eintervalsDistribution[domain]; i < run::mesh->elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return MatrixType::REAL_UNSYMMETRIC;
				}
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void HeatTransferControler::analyticRegularization(size_t domain, bool ortogonalCluster)
{
//	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
//		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
//	}
//
//	if (
//			(_configuration.load_steps_settings.at(_step->step + 1).feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_R &&
//			_configuration.load_steps_settings.at(_step->step + 1).feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_K) ||
//			_configuration.load_steps_settings.at(_step->step + 1).type != LoadStepConfiguration::TYPE::TRANSIENT) {
//
//		if (_configuration.load_steps_settings.at(_step->step + 1).convection.size()) {
//			for (auto it = _configuration.load_steps_settings.at(_step->step + 1).convection.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).convection.end(); ++it) {
//				BoundaryRegionStore *region = _mesh->bregion(it->first);
//				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
//					return;
//				}
//			}
//		}
//
//		if (_configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.size()) {
//			for (auto it = _configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.end(); ++it) {
//				BoundaryRegionStore *region = _mesh->bregion(it->first);
//				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
//					return;
//				}
//			}
//		}
//	}
//
//	if (_configuration.load_steps_settings.at(_step->step + 1).feti.conjugate_projector == FETI_CONJ_PROJECTOR::CONJ_K) {
//		_instance->N1[domain].rows = _instance->K[domain].rows;
//		_instance->N1[domain].cols = 1;
//		_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
//		_instance->N1[domain].type = 'G';
//
//		std::vector<double> diagonal = _instance->K[domain].getDiagonal();
//		_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), diagonal.begin(), diagonal.end());
//
//		_instance->RegMat[domain].rows = _instance->K[domain].rows;
//		_instance->RegMat[domain].cols = _instance->K[domain].cols;
//		_instance->RegMat[domain].nnz  = 1;
//		_instance->RegMat[domain].type = _instance->K[domain].type;
//
//		_instance->RegMat[domain].I_row_indices.push_back(1);
//		_instance->RegMat[domain].J_col_indices.push_back(1);
//		_instance->RegMat[domain].V_values.push_back(_instance->K[domain].getDiagonalMaximum());
//		_instance->RegMat[domain].ConvertToCSR(1);
//	} else {
//		double value;
//		if (ortogonalCluster) {
//			size_t nSum = 0;
//			for (size_t d = 0; d < _instance->domains; d++) {
//				if (_mesh->elements->clusters[d] == _mesh->elements->clusters[domain]) {
//					nSum += _instance->K[d].rows;
//				}
//			}
//			value = 1 / sqrt(nSum);
//		} else {
//			value = 1 / sqrt(_instance->K[domain].rows);
//		}
//
//		_instance->N1[domain].rows = _instance->K[domain].rows;
//		_instance->N1[domain].cols = 1;
//		_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
//		_instance->N1[domain].type = 'G';
//
//		_instance->N1[domain].dense_values.resize(_instance->N1[domain].nnz, value);
//
//		_instance->RegMat[domain].rows = _instance->K[domain].rows;
//		_instance->RegMat[domain].cols = _instance->K[domain].cols;
//		_instance->RegMat[domain].nnz  = 1;
//		_instance->RegMat[domain].type = _instance->K[domain].type;
//
//		_instance->RegMat[domain].I_row_indices.push_back(1);
//		_instance->RegMat[domain].J_col_indices.push_back(1);
//		_instance->RegMat[domain].V_values.push_back(_instance->K[domain].getDiagonalMaximum());
//		_instance->RegMat[domain].ConvertToCSR(1);
//	}
}

MatrixType HeatTransferControler::getMatrixType() const
{
//	if (_step.tangentMatrixCorrection) {
//		return MatrixType::REAL_UNSYMMETRIC;
//	}

	if (_configuration.translation_motions.size()) {
		for (auto it = _configuration.translation_motions.begin(); it != _configuration.translation_motions.end(); ++it) {
			if (run::mesh->eregion(it->first)->eintervals.back().end) {
				return MatrixType::REAL_UNSYMMETRIC;
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

std::vector<double>& HeatTransferControler::getSolutionStore()
{
	return _temperature->data;
}

void HeatTransferControler::dirichletIndices(std::vector<std::vector<eslocal> > &indices)
{
	indices.resize(1); // heat has only one DOF

	for (auto it = _configuration.temperature.regions.begin(); it != _configuration.temperature.regions.end(); ++it) {
		BoundaryRegionStore *region = run::mesh->bregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}

	for (auto it = _configuration.temperature.intersections.begin(); it != _configuration.temperature.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = run::mesh->ibregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}
	_dirichletSize = indices[0].size();
}

void HeatTransferControler::dirichletValues(std::vector<double> &values)
{
	values.resize(_dirichletSize);

	size_t offset = 0;
	for (auto it = _configuration.temperature.regions.begin(); it != _configuration.temperature.regions.end(); ++it) {
		BoundaryRegionStore *region = run::mesh->bregion(it->first);
		it->second.evaluator->evalSelected(
				region->uniqueNodes->datatarray().size(),
				region->uniqueNodes->datatarray().data(),
				3, reinterpret_cast<double*>(run::mesh->nodes->coordinates->datatarray().data()),
				NULL, time::current, values.data() + offset);
		offset += region->uniqueNodes->datatarray().size();
	}

	for (auto it = _configuration.temperature.intersections.begin(); it != _configuration.temperature.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = run::mesh->ibregion(it->first);
		it->second.evaluator->evalSelected(
				region->uniqueNodes->datatarray().size(),
				region->uniqueNodes->datatarray().data(),
				3, reinterpret_cast<double*>(run::mesh->nodes->coordinates->datatarray().data()),
				NULL, time::current, values.data() + offset);
		offset += region->uniqueNodes->datatarray().size();
	}
}


