#include <QDebug>
#include <QApplication>
#include <stdio.h>

#include "mesh/meshwidget.h"
#include "mesh/regionpickerwidget.h"
#include "data/datatype.h"
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

#include "mpi/mpi.h"

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
