
#include "heattransfer2d.controler.h"

#include "../../step.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../basis/matrices/matrixtype.h"
#include "../../../config/ecf/physics/heattransfer.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"
#include "../../../mesh/store/boundaryregionstore.h"

using namespace espreso;

MatrixType HeatTransferControler::getMatrixType() const
{
	if (_step.tangentMatrixCorrection) {
		return MatrixType::REAL_UNSYMMETRIC;
	}

	if (_stepSettings.translation_motions.size()) {
		for (auto it = _stepSettings.translation_motions.begin(); it != _stepSettings.translation_motions.end(); ++it) {
			if (_mesh.eregion(it->first)->eintervals.back().end) {
				return MatrixType::REAL_UNSYMMETRIC;
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

MatrixType HeatTransferControler::getMatrixType(size_t domain) const
{
	if (_step.tangentMatrixCorrection) {
		return MatrixType::REAL_UNSYMMETRIC;
	}

	if (_stepSettings.translation_motions.size()) {
		for (auto it = _stepSettings.translation_motions.begin(); it != _stepSettings.translation_motions.end(); ++it) {
			ElementsRegionStore *region = _mesh.eregion(it->first);
			for (eslocal i = _mesh.elements->eintervalsDistribution[domain]; i < _mesh.elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return MatrixType::REAL_UNSYMMETRIC;
				}
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

	for (auto it = _stepSettings.temperature.regions.begin(); it != _stepSettings.temperature.regions.end(); ++it) {
		BoundaryRegionStore *region = _mesh.bregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}

	for (auto it = _stepSettings.temperature.intersections.begin(); it != _stepSettings.temperature.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = _mesh.ibregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}
	_dirichletSize = indices[0].size();
}

void HeatTransferControler::dirichletValues(std::vector<double> &values)
{
	values.resize(_dirichletSize);

	size_t offset = 0;
	for (auto it = _stepSettings.temperature.regions.begin(); it != _stepSettings.temperature.regions.end(); ++it) {
		BoundaryRegionStore *region = _mesh.bregion(it->first);
		it->second.evaluator->evalSelected(
				region->uniqueNodes->datatarray().size(),
				region->uniqueNodes->datatarray().data(),
				3, reinterpret_cast<double*>(_mesh.nodes->coordinates->datatarray().data()),
				NULL, _step.currentTime, values.data() + offset);
		offset += region->uniqueNodes->datatarray().size();
	}

	for (auto it = _stepSettings.temperature.intersections.begin(); it != _stepSettings.temperature.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = _mesh.ibregion(it->first);
		it->second.evaluator->evalSelected(
				region->uniqueNodes->datatarray().size(),
				region->uniqueNodes->datatarray().data(),
				3, reinterpret_cast<double*>(_mesh.nodes->coordinates->datatarray().data()),
				NULL, _step.currentTime, values.data() + offset);
		offset += region->uniqueNodes->datatarray().size();
	}
}


