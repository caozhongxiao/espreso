#include "declarationswidget.h"
#include "ui_declarationswidget.h"

#include "../models/treemodel.h"
#include "declarations/variabledialog.h"
#include "../data/variable.h"
#include "declarations/materialdialog.h"

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
    this->treeNodeVars = new QStandardItem(tr("Variables"));
    this->treeNodeCS = new QStandardItem(tr("Coordinate Systems"));
    this->treeNodeMats = new QStandardItem(tr("Materials"));
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
    treeMenu.addAction(this->editItem);
    treeMenu.addAction(this->delItem);
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

    QModelIndex clicked = indexList.at(0);
    QModelIndex parent = clicked.parent();

    if (!parent.isValid())
    {
        parent = clicked;
    }

    if (parent.row() == 0)
    {
        VariableDialog* dialog = new VariableDialog(this->varDict, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            Variable v = dialog->data();
            this->variables.append(v);
            this->varDict.insert(v.name(), v);
            QStandardItem* item = new QStandardItem(v.toString());
            treeNodeVars->appendRow(item);
        }
    }
    else if (parent.row() == 1)
    {

    }
    else if (parent.row() == 2)
    {
        MaterialDialog* dialog = new MaterialDialog(this->matDict, this->varDict, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            Material m = dialog->data();
            this->materials.append(m);
            this->matDict.insert(m.name(), m);
            QStandardItem* item = new QStandardItem(m.toString());
            treeNodeMats->appendRow(item);
        }
    }
    else {
        qWarning("%s", QString(tr("Unknown item in declarations!")).toStdString().c_str());
    }
}

void DeclarationsWidget::treeEditItem()
{
    QModelIndexList indexList = ui->DeclarationTree->selectionModel()->selectedIndexes();
    if (!indexList.count())
        return;

    QModelIndex clicked = indexList.at(0);

    this->createEditDialog(clicked);
}

void DeclarationsWidget::treeDelItem()
{
    QModelIndexList indexList = ui->DeclarationTree->selectionModel()->selectedIndexes();
    if (!indexList.count())
        return;

    QModelIndex clicked = indexList.at(0);
    QModelIndex parent = clicked.parent();

    if (!parent.isValid())
        return;

    if (parent.row() == 0)
    {
        Variable v = this->variables.at(clicked.row());
        this->varDict.remove(v.name());
        this->variables.remove(clicked.row());
        treeNodeVars->removeRow(clicked.row());
    }
    else if (parent.row() == 1)
    {

    }
    else if (parent.row() == 2)
    {

    }
    else {
        qWarning("%s", QString(tr("Unknown item in declarations!")).toStdString().c_str());
    }
}

void DeclarationsWidget::createEditDialog(const QModelIndex& item)
{
    QModelIndex parent = item.parent();

    if (!parent.isValid())
        return;

    if (parent.row() == 0)
    {
        VariableDialog* dialog = new VariableDialog(this->variables.at(item.row()),
                                                    this->varDict, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            Variable v = dialog->data();
            this->variables[item.row()] = v;
            this->varDict[v.name()] = v;
            QStandardItem* editted = new QStandardItem(v.toString());
            treeNodeVars->insertRow(item.row() + 1, editted);
            treeNodeVars->removeRow(item.row());
        }
    }
    else if (parent.row() == 1)
    {

    }
    else if (parent.row() == 2)
    {
        MaterialDialog* dialog = new MaterialDialog(this->materials.at(item.row()),
                                                    this->matDict, this->varDict, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            Material m = dialog->data();
            this->materials[item.row()] = m;
            this->matDict[m.name()] = m;
            QStandardItem* editted = new QStandardItem(m.toString());
            treeNodeMats->insertRow(item.row() + 1, editted);
            treeNodeMats->removeRow(item.row());
        }
    }
    else {
        qWarning("%s", QString(tr("Unknown item in declarations!")).toStdString().c_str());
    }
}

void DeclarationsWidget::on_DeclarationTree_doubleClicked(const QModelIndex &index)
{
    this->createEditDialog(index);
}
