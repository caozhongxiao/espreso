#ifndef MESHWIDGET_H
#define MESHWIDGET_H

#include <QtGui>
#include <QOpenGLWidget>

#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/region.h"
#include "../mesh/elements/plane/planeelement.h"

namespace espreso {

    struct MeshPoint3D
    {
        eslocal id;
        QVector3D vector;
    };

    struct MeshTriangle
    {
        QVector<MeshPoint3D> points;
    };

    struct MeshFace
    {
        QVector<MeshTriangle> triangles;
    };

    struct MeshElement
    {
        QVector<MeshFace> faces;
    };

    struct MeshRegion
    {
        QVector<float> points;
        QVector3D color;
        bool isActive = true;
    };

    class MeshWidget : public QOpenGLWidget, protected QOpenGLFunctions
    {
        Q_OBJECT

    public:
        MeshWidget(QWidget* parent = 0);
        MeshWidget(Mesh* mesh, QWidget* parent = 0);
        ~MeshWidget();

        static void initOGL();

        QList<QString> regions();
        void changeRegionState(const QString&);

    signals:
        void regionClicked(const QString& region);

    protected:
        virtual void initializeGL() override;
        virtual void paintGL() override;
        virtual void resizeGL(int w, int h) override;

        virtual void mouseMoveEvent(QMouseEvent *event) override;
        virtual void wheelEvent(QWheelEvent *event) override;
        virtual void mousePressEvent(QMouseEvent *event) override;

    private:
        Mesh* m_mesh;

        float* m_triangles;
        QVector<MeshElement> m_elements;
        QVector<QVector<float> > m_rawelements;
        QMap<QString, MeshRegion> m_regions;
        void computeMesh();

        QOpenGLShaderProgram* m_basicProgram;
        GLuint m_basicProgram_position;
        GLuint m_basicProgram_normal;
        QVector3D m_basicProgram_objectColor;
        QVector3D m_basicProgram_lightColor;
        QVector3D m_basicProgram_lightPos;

        float m_lastX;
        float m_lastY;
        float m_viewRotX = 0;
        float m_viewRotY = 0;
        float m_mouse = true;

        float m_fov = 45.0f;

        QOpenGLShaderProgram* m_clickProgram;
        GLuint m_clickProgram_position;

        const char* m_basicVS =
                "#version 330 core \n"
                "attribute highp vec3 aPos;\n"
                "attribute highp vec3 aNormal;\n"
                "uniform mat4 view;\n"
                "uniform mat4 projection;\n"
                "out vec3 Normal;"
                "out vec3 FragPos;"
                "void main()\n"
                "{\n"
                "gl_Position = projection * view * vec4(aPos, 1.0f);\n"
                "FragPos = aPos;\n"
                "Normal = aNormal;\n"
                "}\n";

        const char* m_basicFS =
                "#version 330 core \n"
                "in vec3 Normal;\n"
                "in vec3 FragPos;\n"
                "uniform vec3 objectColor;\n"
                "uniform vec3 lightColor;\n"
                "uniform vec3 lightPos;\n"
                "out vec4 FragColor;\n"
                "void main()\n"
                "{\n"
                "float ambientStrength = 0.5;\n"
                "vec3 ambient = ambientStrength * lightColor;\n"
                "vec3 norm = normalize(Normal);\n"
                "vec3 lightDir = normalize(lightPos - FragPos);\n"
                "float diff = max(dot(norm, lightDir), 0.0f);\n"
                "vec3 diffuse = diff * lightColor;\n"
                "vec3 result = (ambient + diffuse) * objectColor;\n"
                "FragColor = vec4(result, 1.0f);\n"
                "}\n";

        const char* m_clickVS =
                "#version 330 core\n"
                "uniform mat4 projection;\n"
                "uniform mat4 view;\n"
                "in vec4 position;\n"
                "void main()\n"
                "{\n"
                "gl_Position = projection * view * position;"
                "}\n";
        const char* m_clickFS =
                "#version 330\n"
                "uniform vec3 code;\n"
                "out vec4 outputF;\n"
                "void main()\n"
                "{\n"
                "outputF = vec4(code.x / 255.0f, code.y / 255.0f, code.z / 255.0f, 0);"
                "}\n";
    };

}

#endif // MESHWIDGET_H
