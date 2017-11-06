#ifndef WORKSPACEWINDOW_H
#define WORKSPACEWINDOW_H

#include <QMainWindow>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"
#include "../../input/loader.h"

#include "../declarations/declarationswidget.h"
#include "../mesh/meshwidget.h"
#include "../mesh/regionpickerwidget.h"
#include "workflowwidget.h"

#include "../parallel/mpimanager.h"

namespace espreso
{

namespace Ui {
class WorkspaceWindow;
}

class WorkspaceWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit WorkspaceWindow(MpiManager* manager, QWidget *parent = 0);
    ~WorkspaceWindow();

    void init();

private slots:
    void onFileOpened(const QString& filename);
    void onPhysicsChanged(ECFObject* physics);

private:
    Ui::WorkspaceWindow *ui;

    MpiManager* m_manager;

    ECFConfiguration* m_ecf;
    Mesh* m_mesh;

    WorkflowWidget* m_workflow;

    DeclarationsWidget* m_declarations = nullptr;
    MeshWidget* m_mesh3D = nullptr;
    RegionPickerWidget* m_regions = nullptr;

    void initUi();
    void initPanels();
};

}

#endif // WORKSPACEWINDOW_H
