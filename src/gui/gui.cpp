#include <QDebug>
#include <QApplication>

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

    ECFConfiguration ecf;

    MpiManager mpim;

    WorkspaceWindow ww;

    if (environment->MPIrank == 0) ww.show();

    if (environment->MPIrank != 0)
    {
        while (true)
        {

        }
    }

//    ECFConfiguration ecf(&argc, &argv);
//    Mesh mesh;
//    input::Loader::load(ecf, mesh, environment->MPIrank, environment->MPIsize);

//    if (environment->MPIrank == 0) MeshWidget::initOGL();
//    MeshWidget w(&mesh);
//    RegionPickerWidget rpw(&w);

//    if (environment->MPIrank == 0)
//    {
//        w.show();
//        rpw.show();
//    }


//    for (size_t e = 0; e < mesh.elements().size(); e++) {
//        mesh.elements()[e]->fillFaces();
//    }

//    QMap<QString, QVector<float> > regions;

//    for (size_t e = 0; e < mesh.elements().size(); e++) {

//        for (size_t f = 0; f < mesh.elements()[e]->faces(); f++) {
//            QVector<QString> regionNames;
//            regionNames << QLatin1String("#global");
//            if (mesh.elements()[e]->face(f)->regions().size())
//            {
//                regionNames.clear();

//                for (size_t r = 0; r < mesh.elements()[e]->face(f)->regions().size(); r++)
//                {
//                    QString regionName = QString::fromStdString(mesh.elements()[e]->face(f)->regions()[0]->name);
//                    regionNames << regionName;

//                    if (!regions.contains(regionName))
//                    {
//                        regions.insert(regionName, QVector<float>());
//                    }
//                }
//            }

//            std::vector<std::vector<eslocal> > triangles = dynamic_cast<PlaneElement*>(mesh.elements()[e]->face(f))->triangularize();

//            for (size_t t = 0; t < triangles.size(); t++) {

//                for (size_t n = 0; n < triangles[t].size(); n++) {
//                    foreach (QString rn, regionNames)
//                    {
//                        regions[rn].push_back(mesh.coordinates()[triangles[t][n]].x);
//                        regions[rn].push_back(mesh.coordinates()[triangles[t][n]].y);
//                        regions[rn].push_back(mesh.coordinates()[triangles[t][n]].z);
//                        regions[rn].push_back(0.0f);
//                        regions[rn].push_back(1.0f);
//                        regions[rn].push_back(0.0f);
//                    }
//                }

//            }
//        }
//    }

//    QMapIterator<QString, QVector<float> > it(regions);

//    QMap<QString, QVector<float> > gatheredRegions;

//    while (it.hasNext())
//    {
//        it.next();

//        int num = it.value().size();
//        int nums[environment->MPIsize];
//        MPI_Gather(&num, 1, MPI_INT, nums, 1, MPI_INT, 0, environment->MPICommunicator);

//        int numsum = 0;
//        int displs[environment->MPIsize];
//        QVector<float> coordinates;
//        if (environment->MPIrank == 0)
//        {
//            for (int i = 0; i < environment->MPIsize; i++)
//            {
//                displs[i] = numsum;
//                numsum += nums[i];
//            }
//            coordinates.resize(numsum);
//        }

//        MPI_Gatherv(it.value().data(), num, MPI_FLOAT, coordinates.data(), nums, displs, MPI_FLOAT, 0, environment->MPICommunicator);

//        if (environment->MPIrank == 0)
//        {
//            gatheredRegions.insert(it.key(), coordinates);
//        }
//    }

    if (environment->MPIrank == 0) a.exec();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
