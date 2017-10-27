#include "workflowwidget.h"
#include "ui_workflowwidget.h"

#include <QFileDialog>
#include <QLabel>
#include <QComboBox>
#include <QDebug>
#include <QScrollArea>

#include "loadstepwidget.h"

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

    this->createMaterialsTab();

    this->createLoadstepsTabs();

    this->m_loadsteps_fst_tab_index = ui->workflow->count();
}

void WorkflowWidget::createPhysicsTab()
{
    PhysicsWidget* pw = new PhysicsWidget(this->m_ecf, this);
    pw->init();

    this->m_phyDetail = pw;
    this->m_loadsteps = QString::fromStdString(
                pw->activePhysics()
                        ->getParameter("load_steps")
                        ->getValue()
                ).toInt();
    connect(pw, &PhysicsWidget::loadstepsChanged, this, &WorkflowWidget::onLoadstepsChange);

    ui->workflow->addTab(pw, QLatin1String("Physics"));
    ui->workflow->setCurrentIndex(ui->workflow->count() - 1);
}

void WorkflowWidget::createMaterialsTab()
{

}

void WorkflowWidget::createLoadstepsTabs()
{
    for (int i = 0; i < m_loadsteps; i++)
    {
        LoadstepWidget* lsw = new LoadstepWidget(i + 1, m_phyDetail->activePhysics(), this);
        lsw->init();
        ui->workflow->addTab(lsw, tr("Loadstep %1").arg(i + 1));
    }
}

void WorkflowWidget::onLoadstepsChange(int loadsteps)
{
    int delta = this->m_loadsteps - loadsteps;
    if (delta > 0)
    {
        //DELETE LOADSTEPS

        ui->workflow->removeTab(this->m_loadsteps + this->m_loadsteps_fst_tab_index - 2);

        this->m_loadsteps--;
    }
    else if (delta < 0)
    {
        //ADD LOADSTEPS

        LoadstepWidget* lsw = new LoadstepWidget(++this->m_loadsteps, m_phyDetail->activePhysics(), this);
        lsw->init();
        ui->workflow->addTab(lsw, tr("Loadstep %1").arg(this->m_loadsteps));
    }
}
