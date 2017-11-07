#ifndef WORKFLOWWIDGET_H
#define WORKFLOWWIDGET_H

#include <QWidget>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/plane/planeelement.h"

#include "physicswidget.h"

namespace espreso
{

namespace Ui {
class WorkflowWidget;
}

class WorkflowWidget : public QWidget
{
    Q_OBJECT

public:
    explicit WorkflowWidget(QWidget *parent = 0);
    ~WorkflowWidget();

    void setData(ECFConfiguration* ecf, Mesh* mesh);
    PhysicsConfiguration* activePhysics(ECFConfiguration* ecf);

signals:
    void fileOpened(const QString& filename);
    void physicsChanged(ECFObject* physics);

private slots:
    void on_btnMesh_clicked();
    void onLoadstepsChange(int loadsteps);
    void onPhysicsChange(ECFObject* physics);

private:
    Ui::WorkflowWidget *ui;

    ECFConfiguration* m_ecf;
    Mesh* m_mesh;

    void createPhysicsTab();
    QWidget* m_physicsTab;
    PhysicsWidget* m_phyDetail = nullptr;

    void createMaterialsTab();

    int m_loadsteps;
    int m_loadsteps_fst_tab_index;
    void createLoadstepsTabs();
};

}

#endif // WORKFLOWWIDGET_H
