#include "forms/modelwidget.h"
#include "forms/workflowwidget.h"
#include "forms/declarationswidget.h"
#include "forms/plot/plot.h"
#include "forms/declarations/material/materialpropertieswidget.h"
#include "forms/declarations/datatypeeditwidget.h"
#include "data/datatype.h"
#include <QDebug>
#include <QApplication>

#include "../config/ecf/environment.h"
#include "../config/ecf/ecf.h"
#include "../config/valueholder.h"
#include "../config/ecf/physics/physics.h"
#include "../basis/expression/expression.h"

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

//    std::string exp{"switch {"
//                    "case x < 0 : x;"
//                    "case x >= 0 : x;"
//                    "default : 0;"
//                    "}"};
//    std::string exp1{"if (x > 0 and x < 4) -1;"};
//    std::vector<std::string> vars;
//    vars.push_back("x");
//    Expression e(exp1, vars);
//    for (int i = 0; i < 5; ++i)
//    {
//        std::vector<double> neg{-i};
//        std::vector<double> pos{i};
//        double negV = e.evaluate(neg);
//        double posV = e.evaluate(pos);
//        std::cout << negV << " " << posV << std::endl;
//    }

//    ECFConfiguration ecf;
//    printECF(ecf, 0);

    MaterialConfiguration mat;
//    for (auto p = mat.thermal_properties.parameters.begin(); p != mat.thermal_properties.parameters.end(); ++p)
//    {
//        std::cout << (*p)->name << std::endl;
//    }
    TensorPropertyWidget tpw(&mat.thermal_properties);
    tpw.show();

//    std::string val = "if (x > 0 and x < 4) -1;";
//    ECFValueHolder<std::string> par(val);
//    par.name = "dens";
//    par.metadata.description.push_back("Density");
//    par.metadata.unit = "kg/m^3";


//    MaterialPropertyTableWidget mptw;
//    mptw.addProperty(par);
//    mptw.show();

//    std::string ifst = "if (x > 0 and x < 4) -1;";
//    ECFValueHolder<std::string> expr(ifst);
//    DataTypeEditWidget w(expr);
//    w.show();
//    std::string switchst = "switch {"
//                           "case x == 1 : 2;"
//                           "case x == 2 : 3;"
//                           "default : 0;"
//                           "}";
//    ECFValueHolder<std::string> expr2(switchst);
//    DataTypeEditWidget w1(expr2);
//    w1.show();
//    std::string expSt = "x^2";
//    ECFValueHolder<std::string> expr3(expSt);
//    DataTypeEditWidget w2(expr3);
//    w2.show();

    a.exec();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
