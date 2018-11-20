
#include "heattransfer2d.controler.h"

#include "../../step.h"

#include "../../../basis/matrices/matrixtype.h"
#include "../../../config/ecf/physics/heattransfer.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"

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




