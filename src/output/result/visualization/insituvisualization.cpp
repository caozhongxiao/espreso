
#include "insituvisualization.h"

#include "../../../wrappers/insituwrappers/insituwrapper.h"

using namespace espreso;

InSituVisualization::InSituVisualization(const Mesh &mesh, const OutputConfiguration &configuration)
: Visualization(mesh, configuration)
{
	_inSitu = new InSituWrapper();
}

InSituVisualization::~InSituVisualization()
{
	delete _inSitu;
}

void InSituVisualization::updateMesh()
{

}

void InSituVisualization::updateSolution(const Step &step)
{
	_inSitu->update();
}




