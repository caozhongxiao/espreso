#include "meshwidget.h"

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

MeshWidget::~MeshWidget()
{
    delete this->m_basicProgram;
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
    this->m_basicProgram_lightPos = QVector3D(1.2f, 1.0f, 2.0f);

    this->m_lastX = width() / 2;
    this->m_lastY = height() / 2;

    glEnable(GL_DEPTH_TEST);
}

void MeshWidget::paintGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
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

    float vertices[] = {
        -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
         0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
         0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
         0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
        -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,

        -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
         0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
        -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,

        -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
        -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
        -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
        -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
        -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
        -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,

         0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
         0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
         0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
         0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
         0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
         0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,

        -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
         0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
         0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
         0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,

        -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
         0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
         0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
        -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
        -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f
    };

    QVector<QVector3D> cubePositions;
    cubePositions << QVector3D(0.0f, 0.0f, 0.0f);

    QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
    vbo.create();
    vbo.bind();
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    vbo.allocate(vertices, sizeof(vertices));

    QOpenGLVertexArrayObject cubeVAO(this);
    cubeVAO.create();
    cubeVAO.bind();

    glVertexAttribPointer(m_basicProgram_position, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(m_basicProgram_normal, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    for (unsigned int i = 0; i < cubePositions.size(); i++)
    {
        QMatrix4x4 model;
        model.translate(cubePositions[i]);
        m_basicProgram->setUniformValue("model", model);

        glDrawArrays(GL_TRIANGLES, 0, 36);
    }

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);

    cubeVAO.release();

    vbo.release();

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
