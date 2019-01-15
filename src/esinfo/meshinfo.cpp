
#include "esinfo/meshinfo.h"
#include "output/result/resultstore.h"

#include <cstddef>

namespace espreso {
namespace info {

Mesh* mesh = NULL;

void storeMesh()
{
	mesh->store->updateMesh();
}

void storeSolution()
{
	mesh->store->updateSolution();
}

}
}






