
#include "insitu.h"

#include "esinfo/ecfinfo.h"

#include "wrappers/catalyst/catalyst.h"

#include <unistd.h>

using namespace espreso;

InSitu::InSitu(const Mesh &mesh)
: Visualization(mesh), _catalyst(NULL)
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
	sleep(info::ecf->output.catalyst_sleep_time);
}




