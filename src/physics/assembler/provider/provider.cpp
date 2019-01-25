
#include "provider.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

Provider::Provider(DataHolder *data, LoadStepConfiguration &configuration)
: _data(data), _configuration(configuration)
{

}

bool Provider::needMatrixVectorProduct()
{
	return _configuration.type == LoadStepConfiguration::TYPE::TRANSIENT;
}

bool Provider::needOriginalStiffnessMatrices()
{
	return _configuration.type == LoadStepConfiguration::TYPE::TRANSIENT;
}

bool Provider::needReactionForces()
{
	return _configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

bool Provider::needOriginalRHS()
{
	return
			_configuration.mode == LoadStepConfiguration::MODE::NONLINEAR &&
			(_configuration.nonlinear_solver.line_search || _configuration.nonlinear_solver.check_second_residual);
}



