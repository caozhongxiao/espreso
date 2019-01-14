
#include "insitu.h"

#include "config/ecf/output.h"

#include "wrappers/catalyst/catalyst.h"

#include <unistd.h>

using namespace espreso;

InSitu::InSitu(const Mesh &mesh, const OutputConfiguration &configuration)
: Visualization(mesh, configuration), _catalyst(NULL)
{

}

InSitu::~InSitu()
{
	if (_catalyst != NULL) {
		delete _catalyst;
	}
}

void InSitu::updateMesh()
{
	_catalyst = new Catalyst();
}

void InSitu::updateSolution()
{
	_catalyst->update();
	sleep(_configuration.catalyst_sleep_time);
}




