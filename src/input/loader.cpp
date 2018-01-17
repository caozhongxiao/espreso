
#include "loader.h"

#include "workbench/workbench.h"

#include "../config/ecf/ecf.h"

#include "../mesh/mesh.h"
#include "../old/input/loader.h"

using namespace espreso;

void Loader::load(const ECFConfiguration &configuration, Mesh &mesh, int MPIrank, int MPIsize)
{
	switch (configuration.input) {
	case INPUT_FORMAT::WORKBENCH:
		WorkbenchLoader::load(configuration, mesh);
		break;
	default:
		input::OldLoader::load(configuration, *mesh.mesh, MPIrank, MPIsize);
		break;
	}

	exit(0);
	mesh.load();
}


