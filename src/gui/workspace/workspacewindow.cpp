#include "workspacewindow.h"
#include "ui_workspacewindow.h"

#include <QFileDialog>

using namespace espreso;

WorkspaceWindow::WorkspaceWindow(MpiManager* manager, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::WorkspaceWindow)
{
    this->m_manager = manager;
}

WorkspaceWindow::~WorkspaceWindow()
{
    delete ui;
}

void WorkspaceWindow::init()
{
    ui->setupUi(this);
    this->initUi();

    this->m_ecf = this->m_manager->ecf();
    this->initPanels();
    this->m_mesh = this->m_manager->mesh();
    this->m_workflow->setData(this->m_ecf, this->m_mesh);
}

void WorkspaceWindow::initUi()
{
    this->m_workflow = new WorkflowWidget(this);
    ui->right->layout()->addWidget(m_workflow);
    connect(m_workflow, &WorkflowWidget::inputChanged, this, &WorkspaceWindow::onInputChanged);
    connect(m_workflow, &WorkflowWidget::physicsChanged, this, &WorkspaceWindow::onPhysicsChanged);
    if (environment->MPIrank == 0) MeshWidget::initOGL();
}

void WorkspaceWindow::onInputChanged()
{
    this->initPanels();
    this->m_mesh = this->m_manager->mesh();
    this->m_workflow->setData(this->m_ecf, this->m_mesh);
}

void WorkspaceWindow::initPanels()
{
    if (this->m_datasets != nullptr)
    {
        ui->left->layout()->removeWidget(m_datasets);
        this->m_datasets = nullptr;
    }

    if (this->m_mesh3D != nullptr)
    {
        ui->center->layout()->removeWidget(m_mesh3D);
        this->m_mesh3D = nullptr;
    }

    if (this->m_regions != nullptr)
    {
        ui->left->layout()->removeWidget(m_regions);
        this->m_regions = nullptr;
    }

    ECFObject* materials = static_cast<ECFObject*>(this->m_workflow->activePhysics(this->m_ecf)->getParameter("materials"));
    this->m_datasets = new DataSetsWidget(materials, tr("Data sets"), this);
    ui->left->layout()->addWidget(m_datasets);

    this->m_mesh3D = new MeshWidget(m_manager, this);
    ui->center->layout()->addWidget(m_mesh3D);
    m_mesh3D->setVisible(true);

    this->m_regions = new RegionPickerWidget(m_mesh3D, this);
    ui->left->layout()->addWidget(m_regions);
}

void WorkspaceWindow::onPhysicsChanged(ECFObject *physics)
{
    ECFObject* materials = static_cast<ECFObject*>(physics->getParameter("materials"));
    this->m_datasets->setMaterials(materials);
}

void espreso::WorkspaceWindow::on_btnOpen_pressed()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Open ESPRESO Configuration File"), ".", tr("ECF (*.ecf)"));
    if (filename.isEmpty()) return;

    this->m_manager->masterOpenECF(filename);

    this->m_ecf = this->m_manager->ecf();
    this->initPanels();
    this->m_mesh = this->m_manager->mesh();
    this->m_workflow->setData(this->m_ecf, this->m_mesh);
}

void espreso::WorkspaceWindow::on_btnSave_pressed()
{

}
