
#include "insituvisualization.h"

#include "../../../config/ecf/output.h"

#include "../../../wrappers/insituwrappers/insituwrapper.h"

#include <unistd.h>

using namespace espreso;

InSituVisualization::InSituVisualization(const Mesh &mesh, const OutputConfiguration &configuration)
: Visualization(mesh, configuration), _inSitu(NULL)
{

}

InSituVisualization::~InSituVisualization()
{
	if (_inSitu != NULL) {
		delete _inSitu;
	}
}

void InSituVisualization::updateMesh()
{
	_inSitu = new InSituWrapper(_mesh);
}

void InSituVisualization::updateSolution(const Step &step)
{
	_inSitu->update(step);
	sleep(_configuration.catalyst_sleep_time);
}




