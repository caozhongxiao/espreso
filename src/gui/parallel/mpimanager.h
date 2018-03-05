#ifndef MPIMANAGER_H
#define MPIMANAGER_H

#include <QtGui>

namespace espreso
{

struct ECFRoot;
class Mesh;


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
    ECFRoot* ecf();

    Mesh* mesh();

    void loop();

    void masterOpenECF(const QString&);
    QMap<QString, QVector<float> >* masterGatherMesh();
    void masterExit();


private:
    ECFRoot* m_ecf;
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
