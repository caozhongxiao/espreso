#include "forms/modelwidget.h"
#include "forms/workflowwidget.h"
#include "forms/datawidget.h"
#include "forms/declarationswidget.h"
#include "forms/plot/plot.h"
#include <QApplication>

#ifdef ESPRESOGUIPROJECTBUILD
#include "../config/ecf/environment.h"
#endif

int main(int argc, char *argv[])
{
#ifdef ESPRESOGUIPROJECTBUILD
    MPI_Init(&argc, &argv);
#endif
    QApplication a(argc, argv);
    ModelWidget model;
    model.show();
    WorkflowWidget workflow;
    workflow.show();
//    DataWidget data;
//    data.show();
    DeclarationsWidget declarations;
    declarations.show();

    Plot plot;
    plot.show();

    a.exec();
#ifdef ESPRESOGUIPROJECTBUILD
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return 0;
}
