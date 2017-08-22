#include "forms/modelwidget.h"
#include "forms/workflowwidget.h"
#include "forms/datawidget.h"
#include "forms/declarationswidget.h"
#include "forms/plot/plot.h"
#include "forms/declarations/material/materialpropertytablewidget.h"
#include "data/datatype.h"
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
//    ModelWidget model;
//    model.show();
//    WorkflowWidget workflow;
//    workflow.show();
//    DeclarationsWidget declarations;
//    declarations.show();

    MaterialPropertyTableWidget mptw;
    mptw.show();
    mptw.addRow("Test", new ExpressionType("1e-12"), "unit", "T");
    QList<QList<QString> > ttData;
    QList<QString> row;
    row << "10" << "20";
    ttData << row;
    mptw.addRow("Test", new TableType(ttData), "unit", "T");
    QList<QList<QString> > ptData;
    QList<QString> r1;
    r1 << "10" << "20" << "x^2";
    ptData << r1;
    mptw.addRow("Test", new PiecewiseFunctionType(ptData), "unit", "t");

//    Plot plot;
//    plot.show();

    a.exec();
#ifdef ESPRESOGUIPROJECTBUILD
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return 0;
}
