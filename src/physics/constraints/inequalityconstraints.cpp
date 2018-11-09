
#include "constraints.h"

#include "../instance.h"
#include "../step.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../mesh/mesh.h"
#include "../../mesh/store/nodestore.h"
#include "../../config/ecf/environment.h"
#include "../../config/expression.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

void Constraints::B1ContactInsert(const Step &step, BoundaryRegionStore *region, const ECFExpressionVector &normals, const ECFExpressionVector &gap)
{
	eslocal LMOffset = _dirichletSize + _gluingSize;
	eslocal d = 0;

	for (int dof = 0; dof < _DOFs; dof++) {
		if (normals.data[dof].isSet()) {
			_instance.B1[d].I_row_indices.insert(_instance.B1[d].I_row_indices.end(), region->nodes->datatarray().size(), 0);
			_instance.B1subdomainsMap[d].insert(_instance.B1subdomainsMap[d].end(), region->nodes->datatarray().size(), 0);
			std::iota(_instance.B1[d].I_row_indices.end() - region->nodes->datatarray().size(), _instance.B1[d].I_row_indices.end(), LMOffset + 1);
			std::iota(_instance.B1subdomainsMap[d].end() - region->nodes->datatarray().size(), _instance.B1subdomainsMap[d].end(), LMOffset);
			for (auto n = region->nodes->datatarray().begin(); n != region->nodes->datatarray().end(); ++n) {
				_instance.B1[d].J_col_indices.push_back(*n * _DOFs + dof + 1);
				_instance.B1[d].V_values.push_back(normals.data[dof].evaluator->evaluate(_mesh.nodes->coordinates->datatarray()[*n], step.timeStep, 0));
				// _instance.B1duplicity[d].push_back(0.5);
				// _instance.B1c[d].push_back(0);
				_instance.B1duplicity[d].push_back(0.5);
				_instance.B1c[d].push_back(gap.data[dof].evaluator->evaluate(_mesh.nodes->coordinates->datatarray()[*n], step.timeStep, 0));
			}
			for (eslocal n = LMOffset; n < LMOffset + (eslocal)region->nodes->datatarray().size(); n++) {
				_instance.B1clustersMap.push_back({ n, environment->MPIrank, (environment->MPIrank + 1) % 2 });
			}
			LMOffset += region->nodes->datatarray().size();
		}
		_instance.B1[d].nnz = _instance.B1[d].I_row_indices.size();
		_instance.B1[d].rows = LMOffset;
		_instance.LB[d].resize(_instance.B1[d].nnz, 0);
	}
	_instance.block[Instance::CONSTRAINT::INEQUALITY_CONSTRAINTS] = _instance.B1[d].rows;
}


