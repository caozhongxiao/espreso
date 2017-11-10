
#include "loader.h"
#include "../newmesh/newmesh.h"


#include "../old/input/loader.h"

using namespace espreso;

void Loader::load(const ECFConfiguration &configuration, NewMesh &mesh, int MPIrank, int MPIsize)
{
	input::OldLoader::load(configuration, mesh.mesh, MPIrank, MPIsize);
	mesh.load();
}


