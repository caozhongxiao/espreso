
#include "loader.h"
#include "../mesh/mesh.h"


#include "../old/input/loader.h"

using namespace espreso;

void Loader::load(const ECFConfiguration &configuration, Mesh &mesh, int MPIrank, int MPIsize)
{
	input::OldLoader::load(configuration, *mesh.mesh, MPIrank, MPIsize);
	mesh.load();
}


