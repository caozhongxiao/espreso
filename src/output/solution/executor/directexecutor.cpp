
#include "directexecutor.h"

#include "../visualization/ensight.h"

#include <iostream>

using namespace espreso;

DirectExecutor::DirectExecutor(const Mesh &mesh, const OutputConfiguration &configuration)
: SolutionStoreExecutor(mesh, configuration), _ensight(NULL)
{
//	switch (configuration.format) {
//	...
//	}
	_ensight = new EnSight("solution", mesh);
}

void DirectExecutor::updateMesh()
{
	_ensight->storeGeometry();
	_ensight->storeFETIData();
}

void DirectExecutor::updateSolution(const Step &step)
{
	_ensight->storeVariables(step);
}

DirectExecutor::~DirectExecutor()
{
	if (_ensight != NULL) {
		delete _ensight;
	}
}



