#ifndef WORKFLOWWIDGET_H
#define WORKFLOWWIDGET_H

#include <QWidget>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

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

    void setECF(ECFConfiguration* ecf);

signals:
    void fileOpened(const QString& filename);

private slots:
    void on_btnMesh_clicked();

private:
    Ui::WorkflowWidget *ui;

    ECFConfiguration* m_ecf;

    void createPhysicsTab();
    QWidget* m_physicsTab;
    PhysicsWidget* m_phyDetail = nullptr;
};

}

#endif // WORKFLOWWIDGET_H
