#include "meshwidget.h"

#include "../../config/ecf/environment.h"
#include "mpi/mpi.h"

using namespace espreso;

void MeshWidget::initOGL()
{
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(format);
}

MeshWidget::MeshWidget(QWidget *parent) : QOpenGLWidget(parent)
{
}

MeshWidget::MeshWidget(Mesh* mesh, QWidget* parent) :
    MeshWidget(parent)
{
    this->m_mesh = mesh;

    this->gatherRegions();

//    this->computeMesh();
}

MeshWidget::~MeshWidget()
{
    delete this->m_basicProgram;
    delete this->m_clickProgram;
}

void MeshWidget::gatherRegions()
{
    for (size_t e = 0; e < m_mesh->elements().size(); e++) {
        m_mesh->elements()[e]->fillFaces();
    }

    QMap<QString, QVector<float> > regions;

    for (size_t e = 0; e < m_mesh->elements().size(); e++) {

        for (size_t f = 0; f < m_mesh->elements()[e]->faces(); f++) {
            QVector<QString> regionNames;
            regionNames << QLatin1String("#global");
            if (m_mesh->elements()[e]->face(f)->regions().size())
            {
                regionNames.clear();

                for (size_t r = 0; r < m_mesh->elements()[e]->face(f)->regions().size(); r++)
                {
                    QString regionName = QString::fromStdString(m_mesh->elements()[e]->face(f)->regions()[0]->name);
                    regionNames << regionName;

                    if (!regions.contains(regionName))
                    {
                        regions.insert(regionName, QVector<float>());
                    }
                }
            }

            std::vector<std::vector<eslocal> > triangles = dynamic_cast<PlaneElement*>(m_mesh->elements()[e]->face(f))->triangularize();

            for (size_t t = 0; t < triangles.size(); t++) {

                for (size_t n = 0; n < triangles[t].size(); n++) {
                    foreach (QString rn, regionNames)
                    {
                        regions[rn].push_back(m_mesh->coordinates()[triangles[t][n]].x);
                        regions[rn].push_back(m_mesh->coordinates()[triangles[t][n]].y);
                        regions[rn].push_back(m_mesh->coordinates()[triangles[t][n]].z);
                        regions[rn].push_back(0.0f);
                        regions[rn].push_back(1.0f);
                        regions[rn].push_back(0.0f);
                    }
                }

            }
        }
    }

    QMapIterator<QString, QVector<float> > it(regions);
//    qint32 colorbits = 10;
    qsrand(environment->MPIrank);

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
            MeshRegion reg;
//            colorbits += 255;
//            float red = (colorbits & 0x00FF0000) >> 16;
//            float green = (colorbits & 0x0000FF00) >> 8;
//            float blue = (colorbits & 0x000000FF);
            float red = (qrand() % 256) / 255.0f;
            float green = qrand() % 256 / 255.0f;
            float blue = qrand() % 256 / 255.0f;
            qInfo() << red << green << blue;
            reg.color = QVector3D(red, green, blue);
            reg.points = coordinates;
            this->m_regions.insert(it.key(), reg);
        }
    }
}

void MeshWidget::computeMesh()
{
    MeshRegion globalRegion;
    globalRegion.color = QVector3D(1.0f, 0.5f, 0.31f);
    this->m_regions.insert(QLatin1String("#global"), globalRegion);

    qint32 colorbits = 0;

    for (size_t e = 0; e < m_mesh->elements().size(); e++) {

        MeshElement element;

        for (size_t f = 0; f < m_mesh->elements()[e]->faces(); f++) {

            QVector<QString> regionNames;
            regionNames << QLatin1String("#global");
            if (m_mesh->elements()[e]->face(f)->regions().size())
            {
                regionNames.clear();

                for (size_t r = 0; r < m_mesh->elements()[e]->face(f)->regions().size(); r++)
                {
                    QString regionName = QString::fromStdString(m_mesh->elements()[e]->face(f)->regions()[0]->name);
                    regionNames << regionName;

                    if (!m_regions.contains(regionName))
                    {
                        MeshRegion mr;
                        colorbits += 10;
                        float red = (colorbits & 0x00FF0000) >> 16;
                        float green = (colorbits & 0x0000FF00) >> 8;
                        float blue = (colorbits & 0x000000FF);
                        mr.color = QVector3D(red, green, blue);
                        m_regions.insert(regionName, mr);
                    }
                }
            }

            MeshFace face;
            std::vector<std::vector<eslocal> > triangles = dynamic_cast<PlaneElement*>(m_mesh->elements()[e]->face(f))->triangularize();

            for (size_t t = 0; t < triangles.size(); t++) {

                MeshTriangle triangle;
                triangle.points << MeshPoint3D() << MeshPoint3D() << MeshPoint3D();

                for (size_t n = 0; n < triangles[t].size(); n++) {
                    foreach (QString rn, regionNames)
                    {
                        m_regions[rn].points
                                    << m_mesh->coordinates()[triangles[t][n]].x
                                    << m_mesh->coordinates()[triangles[t][n]].y
                                    << m_mesh->coordinates()[triangles[t][n]].z;
                        m_regions[rn].points << 0.0f << 1.0f << 0.0f;
                    }
                    triangle.points[n].id = triangles[t][n];
                    triangle.points[n].vector.setX(m_mesh->coordinates()[triangles[t][n]].x);
                    triangle.points[n].vector.setY(m_mesh->coordinates()[triangles[t][n]].y);
                    triangle.points[n].vector.setZ(m_mesh->coordinates()[triangles[t][n]].z);
                }

                face.triangles.append(triangle);
            }

            element.faces.append(face);
        }

        m_elements.append(element);
    }

    foreach (MeshElement e, m_elements) {

        QVector<float> element;

        foreach (MeshFace f, e.faces) {

            foreach (MeshTriangle t, f.triangles) {

                foreach (MeshPoint3D p, t.points) {

                    element << p.vector.x() << p.vector.y() << p.vector.z();
                    element << 0.0f << 1.0f << 0.0f;

                }
            }
        }

        m_rawelements.append(element);
    }
}

void MeshWidget::initializeGL()
{
    this->initializeOpenGLFunctions();

    this->m_basicProgram = new QOpenGLShaderProgram(this);
    m_basicProgram->addShaderFromSourceCode(QOpenGLShader::Vertex, m_basicVS);
    m_basicProgram->addShaderFromSourceCode(QOpenGLShader::Fragment, m_basicFS);
    m_basicProgram->link();

    this->m_basicProgram_position = m_basicProgram->attributeLocation("aPos");
    this->m_basicProgram_normal = m_basicProgram->attributeLocation("aNormal");

    this->m_basicProgram_objectColor = QVector3D(1.0f, 0.5f, 0.31f);
    this->m_basicProgram_lightColor = QVector3D(1.0f, 1.0f, 1.0f);
    this->m_basicProgram_lightPos = QVector3D(1.2f, 1.0f, -2.0f);

    this->m_lastX = width() / 2;
    this->m_lastY = height() / 2;

    this->m_clickProgram = new QOpenGLShaderProgram(this);
    m_clickProgram->addShaderFromSourceCode(QOpenGLShader::Vertex, m_clickVS);
    m_clickProgram->addShaderFromSourceCode(QOpenGLShader::Fragment, m_clickFS);
    m_clickProgram->link();

    this->m_clickProgram_position = m_clickProgram->attributeLocation("position");

    glEnable(GL_DEPTH_TEST);
}

void MeshWidget::paintGL()
{
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_basicProgram->bind();
    m_basicProgram->setUniformValue("objectColor", m_basicProgram_objectColor);
    m_basicProgram->setUniformValue("lightColor", m_basicProgram_lightColor);
    m_basicProgram->setUniformValue("lightPos", m_basicProgram_lightPos);

    QMatrix4x4 view;
    view.translate(0.0f, 0.0f, -3.0f);
    view.rotate(m_viewRotX, 1.0f, 0.0f, 0.0f);
    view.rotate(m_viewRotY, 0.0f, 1.0f, 0.0f);

    m_basicProgram->setUniformValue("view", view);

    QMatrix4x4 projection;
    projection.perspective(m_fov, (float)width() / (float)height(), 0.1f, 100.0f);

    m_basicProgram->setUniformValue("projection", projection);

    foreach (MeshRegion r, m_regions) {

        if (!r.isActive) continue;

        m_basicProgram->setUniformValue("objectColor", r.color);

        float* vertices = &r.points.data()[0];
        int len = r.points.size() / 6;

        QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
        vbo.create();
        vbo.bind();
        vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
        vbo.allocate(vertices, sizeof(float) * r.points.size());

        QOpenGLVertexArrayObject cubeVAO(this);
        cubeVAO.create();
        cubeVAO.bind();

        glVertexAttribPointer(m_basicProgram_position, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(m_basicProgram_normal, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);

        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glDrawArrays(GL_TRIANGLES, 0, len);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        cubeVAO.release();

        vbo.release();
    }

    m_basicProgram->release();
}

void MeshWidget::resizeGL(int w, int h)
{

}

void MeshWidget::mouseMoveEvent(QMouseEvent *event)
{
    float xpos = event->pos().x();
    float ypos = event->pos().y();

    if (m_mouse)
    {
        this->m_lastX = xpos;
        this->m_lastY = ypos;
        m_mouse = false;
    }

    float xoffset = xpos - m_lastX;
    float yoffset = m_lastY - ypos;
    this->m_lastX = xpos;
    this->m_lastY = ypos;

    float sensitivity = 0.35f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    this->m_viewRotX -= yoffset;
    this->m_viewRotY += xoffset;

    if (m_viewRotX > 179.0f)
        m_viewRotX = 179.0f;
    if (m_viewRotX < -179.0f)
        m_viewRotX = -179.0f;

    if (m_viewRotY > 179.0f)
        m_viewRotY = 179.0f;
    if (m_viewRotY < -179.0f)
        m_viewRotY = -179.0f;

    this->update();
}

void MeshWidget::wheelEvent(QWheelEvent *event)
{
    float sign = (event->angleDelta().y() > 0) ? 1.0f : -1.0f;

    if (m_fov >= 1.0f && m_fov <= 45.0f)
        m_fov -= sign * 3;

    if (m_fov <= 1.0f)
        m_fov = 1.0f;

    if (m_fov >= 45.0f)
        m_fov = 45.0f;

    this->update();
}

void MeshWidget::mousePressEvent(QMouseEvent* event)
{
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_clickProgram->bind();

    QMatrix4x4 view;
    view.translate(0.0f, 0.0f, -3.0f);
    view.rotate(m_viewRotX, 1.0f, 0.0f, 0.0f);
    view.rotate(m_viewRotY, 0.0f, 1.0f, 0.0f);

    m_clickProgram->setUniformValue("view", view);

    QMatrix4x4 projection;
    projection.perspective(m_fov, (float)width() / (float)height(), 0.1f, 100.0f);

    m_clickProgram->setUniformValue("projection", projection);

    int color = 1;
    QHash<int, QString> colorCodes;

    QMapIterator<QString, MeshRegion> it(m_regions);

    while (it.hasNext()) {
        it.next();
        if (!it.value().isActive) continue;

        MeshRegion r = it.value();

        float red = (color & 0x00FF0000) >> 16;
        float green = (color & 0x0000FF00) >> 8;
        float blue = (color & 0x000000FF);
        m_clickProgram->setUniformValue("code", QVector3D(red, green, blue));
        colorCodes.insert(color, it.key());
        color++;

        float* vertices = &r.points.data()[0];
        int len = r.points.size() / 6;

        QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
        vbo.create();
        vbo.bind();
        vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
        vbo.allocate(vertices, sizeof(float) * r.points.size());

        QOpenGLVertexArrayObject vao(this);
        vao.create();
        vao.bind();

        glVertexAttribPointer(m_basicProgram_position, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glDrawArrays(GL_TRIANGLES, 0, len);

        glDisableVertexAttribArray(0);

        vao.release();

        vbo.release();
    }

    m_basicProgram->release();

    unsigned char res[4];
    GLint viewport[4];

    int x = event->pos().x();
    int y = event->pos().y();

    glGetIntegerv(GL_VIEWPORT, viewport);
    glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &res);

    int clickedColor = 0;
    clickedColor |= (res[0] << 16);
    clickedColor |= (res[1] << 8);
    clickedColor |= res[2];

    if (colorCodes.contains(clickedColor))
        emit regionClicked(colorCodes[clickedColor]);
}

QList<QString> MeshWidget::regions()
{
    return this->m_regions.keys();
}

void MeshWidget::changeRegionState(const QString& region)
{
    bool _isActive = this->m_regions[region].isActive;
    this->m_regions[region].isActive = !_isActive;

    this->update();
}
