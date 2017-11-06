#ifndef MPIMANAGER_H
#define MPIMANAGER_H

#include <QtGui>

#include "mpi/mpi.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/ecf.h"
#include "../../config/valueholder.h"
#include "../../config/ecf/physics/physics.h"
#include "../../config/ecf/physics/heattransfer.h"
#include "../../basis/expression/expression.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/plane/planeelement.h"

#include "../../input/loader.h"


namespace espreso
{

namespace gui
{
    class MpiOperation
    {
    public:
        static const int EXIT = 0;
        static const int OPEN_ECF = 1;
        static const int GATHER_MESH = 2;
    };
}

class MpiManager
{
public:
    MpiManager(int argc, char* argv[]);
    ~MpiManager();

    bool isECFLoaded() const;
    ECFConfiguration* ecf();

    Mesh* mesh();

    void loop();

    void masterOpenECF(const QString&);
    QMap<QString, QVector<float> >* masterGatherMesh();
    void masterExit();


private:
    ECFConfiguration* m_ecf;
    bool m_ecf_loaded = false;

    Mesh* m_mesh = nullptr;

    void performOperation(int code);

    void slaveOpenECF();
    void slaveGatherMesh();

    void _openECF(const std::string& filename);
    QMap<QString, QVector<float> >* _gatherMesh();
    QVector3D pickColor(float angle);
};

}

#endif // MPIMANAGER_H
