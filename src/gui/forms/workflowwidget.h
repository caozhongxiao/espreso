#ifndef WORKFLOWWIDGET_H
#define WORKFLOWWIDGET_H

#include <QWidget>

namespace Ui {
class WorkflowWidget;
}

class WorkflowWidget : public QWidget
{
    Q_OBJECT

public:
    explicit WorkflowWidget(QWidget *parent = 0);
    ~WorkflowWidget();

private:
    Ui::WorkflowWidget *ui;
};

#endif // WORKFLOWWIDGET_H
