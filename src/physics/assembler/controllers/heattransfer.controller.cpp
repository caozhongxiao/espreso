
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

HeatTransferController::HeatTransferController(int dimension, HeatTransferController *previous, HeatTransferLoadStepConfiguration &configuration)
: Controller(dimension),
  _configuration(configuration),
  _ktemperature(1), _kcoordinate(dimension), _kmotion(dimension), _kheat(1), _kthickness(1),
  _ntemperature(NULL),
  _egradient(NULL), _eflux(NULL), _emotion(NULL), _ethickness(NULL),
  _ephaseChange(NULL), _elatentHeat(NULL)
{
	if (previous) {
		_ntemperature = previous->_ntemperature;
		_egradient = previous->_egradient;
		_eflux = previous->_eflux;
		_emotion = previous->_emotion;
		_ethickness = previous->_ethickness;
		_ephaseChange = previous->_ephaseChange;
		_elatentHeat = previous->_elatentHeat;
	}
}

NodeData* HeatTransferController::solution()
{
	return _ntemperature;
}

void HeatTransferController::dirichletIndices(std::vector<std::vector<esint> > &indices)
{
	indices.resize(1); // heat has only one DOF

	for (auto it = _configuration.temperature.begin(); it != _configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		indices[0].insert(indices[0].end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
	}

	_dirichletSize = indices[0].size();
}

void HeatTransferController::dirichletValues(std::vector<double> &values)
{
	values.resize(_dirichletSize);

	size_t offset = 0;
	for (auto it = _configuration.temperature.begin(); it != _configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		it->second.evaluator->evalSelected(
				region->nodes->datatarray().size(),
				region->nodes->datatarray().data(),
				3, reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()),
				NULL, time::current, values.data() + offset);
		offset += region->nodes->datatarray().size();
	}
}

