
#include "run.h"

#include "../mesh/mesh.h"
#include "../output/result/resultstore.h"

#include <cstddef>

espreso::ECFRoot* espreso::run::ecf = NULL;
espreso::Mesh* espreso::run::mesh = NULL;
espreso::DataHolder* espreso::run::data = NULL;

using namespace espreso;


void run::storeMesh()
{
	run::mesh->store->updateMesh();
}

void run::storeSolution()
{
	run::mesh->store->updateSolution();
}


