
#include "provider.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

Provider::Provider(DataHolder *data, LoadStepSolverConfiguration &configuration)
: _data(data), _configuration(configuration)
{

}

bool Provider::needMatrixVectorProduct()
{
	return _configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT;
}

bool Provider::needOriginalStiffnessMatrices()
{
	return _configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT;
}

bool Provider::needSolverStiffnessMatrices()
{
	return needReactionForces();
}

bool Provider::needReactionForces()
{
	return _configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR;
}

bool Provider::needSolverRHS()
{
	return needReactionForces();
}


