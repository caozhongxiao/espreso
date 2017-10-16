#ifndef WORKFLOWWIDGET_H
#define WORKFLOWWIDGET_H

#include <QWidget>

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

signals:
    void fileOpened(const QString& filename);

private slots:
    void on_btnMesh_clicked();

private:
    Ui::WorkflowWidget *ui;
};

}

#endif // WORKFLOWWIDGET_H
