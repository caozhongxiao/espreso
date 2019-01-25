
#include "esinfo/time.h"
#include "esinfo/meshinfo.h"
#include "heattransfer2d.controller.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "basis/matrices/matrixtype.h"
#include "config/ecf/physics/heattransfer.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"

using namespace espreso;

NodeData* HeatTransferController::solution()
{
	return _temperature;
}

void HeatTransferController::dirichletIndices(std::vector<std::vector<esint> > &indices)
{
	indices.resize(1); // heat has only one DOF

	for (auto it = _configuration.temperature.regions.begin(); it != _configuration.temperature.regions.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}

	for (auto it = _configuration.temperature.intersections.begin(); it != _configuration.temperature.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = info::mesh->ibregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}
	_dirichletSize = indices[0].size();
}

void HeatTransferController::dirichletValues(std::vector<double> &values)
{
	values.resize(_dirichletSize);

	size_t offset = 0;
	for (auto it = _configuration.temperature.regions.begin(); it != _configuration.temperature.regions.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		it->second.evaluator->evalSelected(
				region->uniqueNodes->datatarray().size(),
				region->uniqueNodes->datatarray().data(),
				3, reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()),
				NULL, time::current, values.data() + offset);
		offset += region->uniqueNodes->datatarray().size();
	}

	for (auto it = _configuration.temperature.intersections.begin(); it != _configuration.temperature.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = info::mesh->ibregion(it->first);
		it->second.evaluator->evalSelected(
				region->uniqueNodes->datatarray().size(),
				region->uniqueNodes->datatarray().data(),
				3, reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()),
				NULL, time::current, values.data() + offset);
		offset += region->uniqueNodes->datatarray().size();
	}
}

