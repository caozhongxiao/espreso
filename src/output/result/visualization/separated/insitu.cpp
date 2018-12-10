
#include "insitu.h"

#include "../../../../config/ecf/output.h"

#include "../../../../wrappers/insituwrappers/insituwrapper.h"

#include <unistd.h>

using namespace espreso;

InSitu::InSitu(const Mesh &mesh, const OutputConfiguration &configuration)
: Visualization(mesh, configuration), _inSitu(NULL)
{

}

InSitu::~InSitu()
{
	if (_inSitu != NULL) {
		delete _inSitu;
	}
}

void InSitu::updateMesh()
{
	_inSitu = new InSituWrapper(_mesh);
}

void InSitu::updateSolution()
{
	_inSitu->update();
	sleep(_configuration.catalyst_sleep_time);
}




