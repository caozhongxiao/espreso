#include "workflowwidget.h"
#include "ui_workflowwidget.h"

#include <QFileDialog>

using namespace espreso;

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

void espreso::WorkflowWidget::on_btnMesh_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Open a file"),
                                                    "", tr("Espreso configuration file (*.ecf)"));
    emit fileOpened(filename);
}
