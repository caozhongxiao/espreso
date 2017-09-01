#include "forms/modelwidget.h"
#include "forms/workflowwidget.h"
#include "forms/declarationswidget.h"
#include "forms/plot/plot.h"
#include "forms/declarations/material/materialpropertieswidget.h"
#include "forms/declarations/datatypeeditwidget.h"
#include "data/datatype.h"
#include <QApplication>

#include "../config/ecf/environment.h"
#include "../config/ecf/ecf.h"
#include "../config/valueholder.h"
#include "../config/ecf/physics/physics.h"

using namespace espreso;

template <typename Ttype>
std::ostream& operator<<(std::ostream& os, const std::vector<Ttype> &v)
{
    for (size_t i = 0; i < v.size(); i++) {
        os << v[i] << " ";
    }
    os << "\n";
    return os;
}

void printECF(const ECFObject &object, size_t indent)
{
    auto printindent = [&] () {
        for (size_t j = 0; j < indent; j++) {
            std::cout << " ";
        }
    };

    auto printparameter = [&] (const ECFParameter *parameter) {
        if (parameter == NULL || !parameter->metadata.isallowed()) {
            return;
        }
        printindent();
        std::cout << parameter->name << " ";
        if (parameter->isValue()) {
            std::cout << " = " << parameter->getValue() << ";\n";
            if (parameter->metadata.datatype.front() == ECFDataType::EXPRESSION) {
                printindent();
                std::cout << "EXPRESSION: " << parameter->metadata.variables;
            }
            if (parameter->metadata.tensor != NULL) {
                printindent();
                std::cout << "TENSOR: " << parameter->metadata.tensor->size << "x" << parameter->metadata.tensor->size << "\n";
            }
        } else if (parameter->isObject()) {
            std::cout << "{\n";
            printECF(*dynamic_cast<const ECFObject*>(parameter), indent + 2);
            printindent();
            std::cout << "}\n";
        } else {
            // Separators etc..
            std::cout << "\n";
        }

    };

    for (size_t i = 0; i < object.parameters.size(); i++) {
        printparameter(object.parameters[i]);
    }
};



int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    QApplication a(argc, argv);
//    ModelWidget model;
//    model.show();
//    WorkflowWidget workflow;
//    workflow.show();
//    DeclarationsWidget declarations;
//    declarations.show();

//    TensorProperty p("Thermal Conductivity");
//    QList<TensorPropertyModelItem> model1Items;
//    model1Items << TensorPropertyModelItem("X component", "kg", "KXX", new ExpressionType("10"));
//    TensorPropertyModel model1(1, "Isotropic", model1Items);
//    QList<TensorPropertyModelItem> model2Items;
//    model2Items << TensorPropertyModelItem("X component", "kg", "KXX", new ExpressionType("10"))
//                << TensorPropertyModelItem("Y component", "kg", "KYY", new ExpressionType("10"))
//                << TensorPropertyModelItem("Z component", "kg", "KZZ", new ExpressionType("10"));
//    TensorPropertyModel model2(3, "Diagonal", model2Items);
//    p.appendModel(model1);
//    p.appendModel(model2);

//    ScalarProperty sp("Density", "kJ", "DENS", new ExpressionType("0"));


//    QVector<TensorProperty> tensors;
//    tensors << p;
//    QVector<ScalarProperty> scalars;
//    scalars << sp;
//    MaterialPropertiesWidget w(tensors, scalars);
//    w.show();


//    Plot plot;
//    plot.show();

//    std::string val("x^2");
//    ECFValueHolder<std::string> expr(val);
//    expr.metadata.datatype.push_back(ECFDataType::EXPRESSION);
//    DataTypeEditWidget w(expr);
//    w.show();

    ECFConfiguration ecf(&argc, &argv);
    printECF(ecf, 0);


    a.exec();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
