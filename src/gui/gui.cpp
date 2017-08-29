#include "forms/modelwidget.h"
#include "forms/workflowwidget.h"
#include "forms/datawidget.h"
#include "forms/declarationswidget.h"
#include "forms/plot/plot.h"
#include "forms/declarations/material/materialpropertieswidget.h"
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

    TensorProperty p("Thermal Conductivity");
    QList<TensorPropertyModelItem> model1Items;
    model1Items << TensorPropertyModelItem("X component", "kg", "KXX", new ExpressionType("10"));
    TensorPropertyModel model1(1, "Isotropic", model1Items);
    QList<TensorPropertyModelItem> model2Items;
    model2Items << TensorPropertyModelItem("X component", "kg", "KXX", new ExpressionType("10"))
                << TensorPropertyModelItem("Y component", "kg", "KYY", new ExpressionType("10"))
                << TensorPropertyModelItem("Z component", "kg", "KZZ", new ExpressionType("10"));
    TensorPropertyModel model2(3, "Diagonal", model2Items);
    p.appendModel(model1);
    p.appendModel(model2);

    ScalarProperty sp("Density", "kJ", "DENS", new ExpressionType("0"));


    QVector<TensorProperty> tensors;
    tensors << p;
    QVector<ScalarProperty> scalars;
    scalars << sp;
    MaterialPropertiesWidget w(tensors, scalars);
    w.show();


//    Plot plot;
//    plot.show();

    a.exec();
#ifdef ESPRESOGUIPROJECTBUILD
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return 0;
}
