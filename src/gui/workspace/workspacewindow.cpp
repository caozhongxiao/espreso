#include "workspacewindow.h"
#include "ui_workspacewindow.h"

#include <QFileDialog>

using namespace espreso;

WorkspaceWindow::WorkspaceWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::WorkspaceWindow)
{
    ui->setupUi(this);
    this->initUi();
}

WorkspaceWindow::~WorkspaceWindow()
{
    delete ui;

    if (m_ecf_local) delete this->m_ecf;
}

void WorkspaceWindow::initUi()
{
    this->m_workflow = new WorkflowWidget(this);
    ui->top->layout()->addWidget(m_workflow);
    connect(m_workflow, &WorkflowWidget::fileOpened, this, &WorkspaceWindow::onFileOpened);
    MeshWidget::initOGL();
}

void WorkspaceWindow::onFileOpened(const QString &filename)
{
    this->m_ecf_local = true;
    this->m_ecf = new ECFConfiguration(filename.toStdString());

    this->initPanels();
}

void WorkspaceWindow::initPanels()
{
    ui->left->layout()->addWidget(new DeclarationsWidget(this));
    Mesh mesh;
    input::Loader::load(*m_ecf, mesh, environment->MPIrank, environment->MPIsize);
    MeshWidget* mw = new MeshWidget(&mesh, this);
    ui->center->layout()->addWidget(mw);
    mw->setVisible(true);
    ui->right->layout()->addWidget(new RegionPickerWidget(mw, this));

    qInfo() << environment->MPIrank;
}
