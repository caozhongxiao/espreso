#include "mpimanager.h"

#include <QDebug>

using namespace espreso;

MpiManager::MpiManager(int argc, char* argv[])
{
    this->m_ecf = new ECFConfiguration();
    this->m_ecf_loaded = this->m_ecf->fill(&argc, &argv);
}

bool MpiManager::isECFLoaded() const
{
    return this->m_ecf_loaded;
}

ECFConfiguration* MpiManager::ecf()
{
    return this->m_ecf;
}

void MpiManager::loop()
{
    if (environment->MPIrank == 0) return;

    while (true)
    {
        int opcode = 0;
        MPI_Bcast(&opcode, 1, MPI_INT, 0, environment->MPICommunicator);
        if (opcode == 0) break;

        this->performOperation(opcode);
    }
}

void MpiManager::performOperation(int code)
{
    switch (code)
    {
        case gui::MpiOperation::OPEN_ECF:
            this->slaveOpenECF();
            break;
        case gui::MpiOperation::GATHER_MESH:
            this->slaveGatherMesh();
            break;
        default:
            ;
    }
}

void MpiManager::masterOpenECF(const QString& filename)
{
    int op = gui::MpiOperation::OPEN_ECF;
    MPI_Bcast(&op, 1, MPI_INT, 0, environment->MPICommunicator);

    std::string fname = filename.toStdString();
    std::vector<char> fn(fname.begin(), fname.end());
    fn.push_back(0);
    int len = fname.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, environment->MPICommunicator);
    MPI_Bcast(&fn[0], len, MPI_CHAR, 0, environment->MPICommunicator);

    this->_openECF(filename.toStdString());
}

void MpiManager::slaveOpenECF()
{
    int len = 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, environment->MPICommunicator);
    char fn[len + 1];
    MPI_Bcast(fn, len, MPI_CHAR, 0, environment->MPICommunicator);
    fn[len] = 0;
    this->_openECF(std::string(fn));
}

void MpiManager::_openECF(const std::string& filename)
{
    if (!this->m_ecf_loaded)
    {
        this->m_ecf_loaded = this->m_ecf->fill(filename);
    }
    else
    {
        ECFConfiguration* tmp = this->m_ecf;
        delete tmp;
        this->m_ecf = new ECFConfiguration(filename);
    }
}

QMap<QString, QVector<float> >* MpiManager::masterGatherMesh()
{
    int op = gui::MpiOperation::GATHER_MESH;
    MPI_Bcast(&op, 1, MPI_INT, 0, environment->MPICommunicator);

    return this->_gatherMesh();
}

void MpiManager::slaveGatherMesh()
{
    this->_gatherMesh();
}

QMap<QString, QVector<float> >* MpiManager::_gatherMesh()
{
    Mesh mesh;
    input::Loader::load(*m_ecf, mesh, environment->MPIrank, environment->MPIsize);

    for (size_t e = 0; e < mesh.elements().size(); e++) {
        mesh.elements()[e]->fillFaces();
    }

    QMap<QString, QVector<float> > regions;

    for (size_t e = 0; e < mesh.elements().size(); e++) {

        for (size_t f = 0; f < mesh.elements()[e]->faces(); f++) {
            QVector<QString> regionNames;
            regionNames << QLatin1String("#global");
            if (mesh.elements()[e]->face(f)->regions().size())
            {
                regionNames.clear();

                for (size_t r = 0; r < mesh.elements()[e]->face(f)->regions().size(); r++)
                {
                    QString regionName = QString::fromStdString(mesh.elements()[e]->face(f)->regions()[0]->name);
                    regionNames << regionName;

                    if (!regions.contains(regionName))
                    {
                        regions.insert(regionName, QVector<float>());
                    }
                }
            }

            std::vector<std::vector<eslocal> > triangles = dynamic_cast<PlaneElement*>(mesh.elements()[e]->face(f))->triangularize();

            for (size_t t = 0; t < triangles.size(); t++) {

                for (size_t n = 0; n < triangles[t].size(); n++) {
                    foreach (QString rn, regionNames)
                    {
                        regions[rn].push_back(mesh.coordinates()[triangles[t][n]].x);
                        regions[rn].push_back(mesh.coordinates()[triangles[t][n]].y);
                        regions[rn].push_back(mesh.coordinates()[triangles[t][n]].z);
                        regions[rn].push_back(0.0f);
                        regions[rn].push_back(1.0f);
                        regions[rn].push_back(0.0f);
                    }
                }

            }
        }
    }

    QMapIterator<QString, QVector<float> > it(regions);

    QMap<QString, QVector<float> >* r_regions = (environment->MPIrank == 0) ? new QMap<QString, QVector<float> >() : nullptr;

    while (it.hasNext())
    {
        it.next();

        int num = it.value().size();
        int nums[environment->MPIsize];
        MPI_Gather(&num, 1, MPI_INT, nums, 1, MPI_INT, 0, environment->MPICommunicator);

        int numsum = 0;
        int displs[environment->MPIsize];
        QVector<float> coordinates;
        if (environment->MPIrank == 0)
        {
            for (int i = 0; i < environment->MPIsize; i++)
            {
                displs[i] = numsum;
                numsum += nums[i];
            }
            coordinates.resize(numsum);
        }

        MPI_Gatherv(it.value().data(), num, MPI_FLOAT, coordinates.data(), nums, displs, MPI_FLOAT, 0, environment->MPICommunicator);

        if (environment->MPIrank == 0)
        {
            r_regions->insert(it.key(), coordinates);
        }
    }

    return r_regions;
}

void MpiManager::masterExit()
{
    int code = gui::MpiOperation::EXIT;
    MPI_Bcast(&code, 1, MPI_INT, 0, environment->MPICommunicator);
}
