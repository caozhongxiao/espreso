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

namespace espreso
{

namespace Ui {
class WorkspaceWindow;
}

class WorkspaceWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit WorkspaceWindow(QWidget *parent = 0);
    ~WorkspaceWindow();

private slots:
    void onFileOpened(const QString& filename);

private:
    Ui::WorkspaceWindow *ui;

    bool m_ecf_local = false;
    ECFConfiguration* m_ecf;

    WorkflowWidget* m_workflow;
    void initUi();
    void initPanels();
};

}

#endif // WORKSPACEWINDOW_H
