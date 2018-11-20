
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

TimeStepSolver::TimeStepSolver(const std::string &description, Provider &composer)
: _description(description), _composer(composer)
{

}

std::string TimeStepSolver::description() const
{
	return _description;
}
