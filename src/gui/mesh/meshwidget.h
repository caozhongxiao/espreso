#ifndef MESHWIDGET_H
#define MESHWIDGET_H

#include <QtGui>
#include <QOpenGLWidget>

namespace espreso {

    class MeshWidget : public QOpenGLWidget, protected QOpenGLFunctions
    {
        Q_OBJECT

    public:
        MeshWidget(QWidget* parent = 0);
        ~MeshWidget();

        static void initOGL();

    protected:
        virtual void initializeGL() override;
        virtual void paintGL() override;
        virtual void resizeGL(int w, int h) override;

        virtual void mouseMoveEvent(QMouseEvent *event) override;
        virtual void wheelEvent(QWheelEvent *event) override;

    private:
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

        const char* m_basicVS =
                "#version 330 core \n"
                "attribute highp vec3 aPos;\n"
                "attribute highp vec3 aNormal;\n"
                "uniform mat4 model;\n"
                "uniform mat4 view;\n"
                "uniform mat4 projection;\n"
                "out vec3 Normal;"
                "out vec3 FragPos;"
                "void main()\n"
                "{\n"
                "gl_Position = projection * view * model * vec4(aPos, 1.0f);\n"
                "FragPos = vec3(model * vec4(aPos, 1.0f));\n"
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
    };

}

#endif // MESHWIDGET_H
