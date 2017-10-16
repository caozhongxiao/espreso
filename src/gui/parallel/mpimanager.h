#ifndef MPIMANAGER_H
#define MPIMANAGER_H

#include "mpi/mpi.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/ecf.h"
#include "../../config/valueholder.h"
#include "../../config/ecf/physics/physics.h"
#include "../../config/ecf/physics/heattransfer.h"
#include "../../basis/expression/expression.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/elements/plane/planeelement.h"

#include "../../input/loader.h"


namespace espreso
{

class MpiManager
{
public:
    MpiManager();
};

}

#endif // MPIMANAGER_H
