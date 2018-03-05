#include <QDebug>
#include <QApplication>
#include <stdio.h>

#include "mesh/meshwidget.h"
#include "mesh/regionpickerwidget.h"
#include "parallel/mpimanager.h"

#include "workspace/workspacewindow.h"

#include "../config/ecf/environment.h"
#include "../config/valueholder.h"
#include "../config/ecf/physics/physics.h"
#include "../config/ecf/physics/heattransfer.h"
#include "../basis/expression/expression.h"
#include "../config/ecf/root.h"

#include "../mesh/mesh.h"

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

        // CSS
        QCoreApplication::setAttribute(Qt::AA_UseStyleSheetPropagationInWidgetStyles, true);
        QFile css(":/stylesheets/stylesheets/general.css");
        css.open(QFile::ReadOnly);
        QTextStream stream(&css);
        a.setStyleSheet(stream.readAll());
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
