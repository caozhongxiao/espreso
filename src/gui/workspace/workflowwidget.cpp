#include "workflowwidget.h"
#include "ui_workflowwidget.h"

#include <QFileDialog>
#include <QLabel>
#include <QComboBox>
#include <QDebug>
#include <QScrollArea>

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

void WorkflowWidget::on_btnMesh_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Open a file"),
                                                    "", tr("Espreso configuration file (*.ecf)"));
    if (filename.size() == 0)
        return;

    emit fileOpened(filename);
}

void WorkflowWidget::setECF(ECFConfiguration *ecf)
{
    int tabs = ui->workflow->count();
    for (int i = 1; i < tabs; i++)
    {
        ui->workflow->removeTab(i);
    }

    this->m_ecf = ecf;
    this->m_physicsTab = nullptr;

    this->createPhysicsTab();
}

void WorkflowWidget::createPhysicsTab()
{
    PhysicsWidget* pw = new PhysicsWidget(this->m_ecf, this);
    pw->init();

    ui->workflow->addTab(pw, QLatin1String("Physics"));
    ui->workflow->setCurrentIndex(ui->workflow->count() - 1);
}

