#include <QDebug>
#include <QApplication>
#include <stdio.h>

#include "mesh/meshwidget.h"
#include "mesh/regionpickerwidget.h"
#include "parallel/mpimanager.h"

#include "workspace/workspacewindow.h"

#include "../config/ecf/environment.h"
#include "../config/ecf/ecf.h"
#include "../config/valueholder.h"
#include "../config/ecf/physics/physics.h"
#include "../config/ecf/physics/heattransfer.h"
#include "../basis/expression/expression.h"

#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/elements/plane/planeelement.h"

#include "../input/loader.h"

#include "mpi.h"

using namespace espreso;


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    QApplication a(argc, argv);

    MpiManager mpim(argc, argv);

    WorkspaceWindow ww(&mpim);

    if (environment->MPIrank == 0)
    {
        ww.init();
        ww.show();
    }

    mpim.loop();

    if (environment->MPIrank == 0)
    {
        a.exec();
        mpim.masterExit();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
