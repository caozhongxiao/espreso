#include "declarationswidget.h"
#include "ui_declarationswidget.h"

#include "../models/treemodel.h"
#include "declarations/variabledialog.h"
#include "../data/variable.h"

#include <QStandardItemModel>

DeclarationsWidget::DeclarationsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DeclarationsWidget)
{
    ui->setupUi(this);

//    ui->DeclarationTree->setModel(new TreeModel(tr("Declarations")));
    this->setupTree();
    this->createActions();
}

DeclarationsWidget::~DeclarationsWidget()
{
    delete ui;
}

void DeclarationsWidget::setupTree()
{
    this->treeModel = new QStandardItemModel();
    QStandardItem* parent = this->treeModel->invisibleRootItem();
    this->treeNodeVars = new QStandardItem("Variables");
    this->treeNodeCS = new QStandardItem("Coordinate Systems");
    this->treeNodeMats = new QStandardItem("Materials");
    parent->appendRow(treeNodeVars);
    parent->appendRow(treeNodeCS);
    parent->appendRow(treeNodeMats);
    ui->DeclarationTree->setModel(this->treeModel);
    ui->DeclarationTree->setEditTriggers(QAbstractItemView::NoEditTriggers);
}


void DeclarationsWidget::on_DeclarationTree_customContextMenuRequested(const QPoint &pos)
{
    QMenu treeMenu(this);
    treeMenu.addAction(this->newItem);
    treeMenu.exec(ui->DeclarationTree->mapToGlobal(pos));
}

void DeclarationsWidget::createActions()
{
    this->newItem = new QAction(tr("&New"), this);
    connect(this->newItem, &QAction::triggered, this, &DeclarationsWidget::treeNewItem);

    this->editItem = new QAction(tr("&Edit"), this);
    connect(this->editItem, &QAction::triggered, this, &DeclarationsWidget::treeEditItem);

    this->delItem = new QAction(tr("&Delete"), this);
    connect(this->delItem, &QAction::triggered, this, &DeclarationsWidget::treeDelItem);
}

void DeclarationsWidget::treeNewItem()
{
    QModelIndexList indexList = ui->DeclarationTree->selectionModel()->selectedIndexes();
    if (!indexList.count())
        return;

    VariableDialog* dialog = new VariableDialog(this);
    if (dialog->exec() == QDialog::Accepted)
    {
        QStandardItem* item = new QStandardItem("Test");
        treeNodeVars->appendRow(item);
    }
}

void DeclarationsWidget::treeEditItem()
{

}

void DeclarationsWidget::treeDelItem()
{

}
