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
#include "inputwidget.h"

namespace espreso
{

namespace Ui {
class WorkflowWidget;
}

class WorkflowWidget : public QWidget, public ISavableObject, public IValidatableObject
{
    Q_OBJECT

public:
    explicit WorkflowWidget(QWidget *parent = 0);
    ~WorkflowWidget();

    void setData(ECFConfiguration* ecf, Mesh* mesh);
    PhysicsConfiguration* activePhysics(ECFConfiguration* ecf);
    ECFObject* input();

    virtual void save() override;
    virtual bool isValid() override;
    virtual QString errorMessage() override;

signals:
    void inputChanged();
    void physicsChanged(ECFObject* physics);

private slots:
    void onLoadstepsChange(int loadsteps);
    void onPhysicsChange(ECFObject* physics);
    void onInputChange(int index);

    void on_btnLoad_pressed();

private:
    Ui::WorkflowWidget *ui;

    ECFConfiguration* m_ecf;
    Mesh* m_mesh;

    bool m_inputBox_filled = false;
    void createInput();
    InputWidget* m_inputWidget = nullptr;

    void createPhysicsTab();
    QWidget* m_physicsTab;
    PhysicsWidget* m_phyDetail = nullptr;

    void createMaterialsTab();

    int m_loadsteps;
    int m_loadsteps_fst_tab_index;
    void createLoadstepsTabs();

    void createOutputTab();
};

}

#endif // WORKFLOWWIDGET_H
