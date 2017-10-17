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

    if (this->m_manager->isECFLoaded())
    {
        this->m_ecf = this->m_manager->ecf();
        this->initPanels();
    }
}

void WorkspaceWindow::initUi()
{
    this->m_workflow = new WorkflowWidget(this);
    ui->top->layout()->addWidget(m_workflow);
    connect(m_workflow, &WorkflowWidget::fileOpened, this, &WorkspaceWindow::onFileOpened);
    if (environment->MPIrank == 0) MeshWidget::initOGL();
}

void WorkspaceWindow::onFileOpened(const QString &filename)
{
    this->m_manager->masterOpenECF(filename);

    this->m_ecf = this->m_manager->ecf();

    this->initPanels();
}

void WorkspaceWindow::initPanels()
{
    if (this->m_declarations != nullptr)
    {
        ui->left->layout()->removeWidget(m_declarations);
        this->m_declarations = nullptr;
    }

    if (this->m_mesh != nullptr)
    {
        ui->center->layout()->removeWidget(m_mesh);
        this->m_mesh = nullptr;
    }

    if (this->m_regions != nullptr)
    {
        ui->right->layout()->removeWidget(m_regions);
        this->m_regions = nullptr;
    }

    this->m_declarations = new DeclarationsWidget(this);
    ui->left->layout()->addWidget(m_declarations);

    this->m_mesh = new MeshWidget(m_manager, this);
    ui->center->layout()->addWidget(m_mesh);
    m_mesh->setVisible(true);

    this->m_regions = new RegionPickerWidget(m_mesh, this);
    ui->right->layout()->addWidget(m_regions);
}
