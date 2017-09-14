#include "workflowwidget.h"
#include "ui_workflowwidget.h"

WorkflowWidget::WorkflowWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WorkflowWidget)
{
    ui->setupUi(this);
}

WorkflowWidget::~WorkflowWidget()
{
    delete ui;
}
