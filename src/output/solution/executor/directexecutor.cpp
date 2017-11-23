
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
}

void DirectExecutor::updateSolution()
{
	_ensight->storeVariables();
}

DirectExecutor::~DirectExecutor()
{
	if (_ensight != NULL) {
		delete _ensight;
	}
}



