#include "declarationswidget.h"
#include "ui_declarationswidget.h"

#include "declarations/variabledialog.h"
#include "../data/variable.h"
#include "declarations/material/materialdialog.h"
#include "../config/ecf/material/material.h"

#include <QStandardItemModel>

using namespace espreso;

DeclarationsWidget::DeclarationsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DeclarationsWidget)
{
    ui->setupUi(this);

    this->setupTree();
    this->createActions();
}

DeclarationsWidget::~DeclarationsWidget()
{
    delete ui;
}

void DeclarationsWidget::setupTree()
{
    this->m_treeModel = new QStandardItemModel();
    QStandardItem* parent = this->m_treeModel->invisibleRootItem();
    //this->m_treeNodeVars = new QStandardItem(tr("Variables"));
    //this->m_treeNodeCS = new QStandardItem(tr("Coordinate Systems"));
    this->m_treeNodeMats = new QStandardItem(tr("Materials"));

//    parent->appendRow(m_treeNodeVars);
//    this->m_varRow = m_treeModel->rowCount() - 1;
//    parent->appendRow(m_treeNodeCS);
//    this->m_csRow = m_treeModel->rowCount() - 1;
    parent->appendRow(m_treeNodeMats);
    this->m_matRow = m_treeModel->rowCount() - 1;

    ui->DeclarationTree->setModel(this->m_treeModel);
    ui->DeclarationTree->setEditTriggers(QAbstractItemView::NoEditTriggers);
}


void DeclarationsWidget::on_DeclarationTree_customContextMenuRequested(const QPoint &pos)
{
    QMenu treeMenu(this);
    treeMenu.addAction(this->m_newItem);
    treeMenu.addAction(this->m_editItem);
    treeMenu.addAction(this->m_delItem);
    treeMenu.exec(ui->DeclarationTree->mapToGlobal(pos));
}

void DeclarationsWidget::createActions()
{
    this->m_newItem = new QAction(tr("&New"), this);
    connect(this->m_newItem, &QAction::triggered, this, &DeclarationsWidget::treeNewItem);

    this->m_editItem = new QAction(tr("&Edit"), this);
    connect(this->m_editItem, &QAction::triggered, this, &DeclarationsWidget::treeEditItem);

    this->m_delItem = new QAction(tr("&Delete"), this);
    connect(this->m_delItem, &QAction::triggered, this, &DeclarationsWidget::treeDelItem);
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

    if (parent.row() == m_varRow)
    {
        VariableDialog* dialog = new VariableDialog(this->m_varDict, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            Variable v = dialog->data();
            this->m_variables.append(v);
            this->m_varDict.insert(v.name(), v);
            QStandardItem* item = new QStandardItem(v.toString());
            m_treeNodeVars->appendRow(item);
        }
    }
    else if (parent.row() == m_csRow)
    {

    }
    else if (parent.row() == m_matRow)
    {
        MaterialConfiguration material;
        MaterialDialog* dialog = new MaterialDialog(&material, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            QStandardItem* item = new QStandardItem(
                        QString::fromStdString(material.getParameter("name")->getValue())
                        );
            m_treeNodeMats->appendRow(item);
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

    if (parent.row() == m_varRow)
    {
        Variable v = this->m_variables.at(clicked.row());
        this->m_varDict.remove(v.name());
        this->m_variables.remove(clicked.row());
        m_treeNodeVars->removeRow(clicked.row());
    }
    else if (parent.row() == m_csRow)
    {

    }
    else if (parent.row() == m_matRow)
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

    if (parent.row() == m_varRow)
    {
        VariableDialog* dialog = new VariableDialog(this->m_variables.at(item.row()),
                                                    this->m_varDict, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            Variable v = dialog->data();
            this->m_variables[item.row()] = v;
            this->m_varDict[v.name()] = v;
            QStandardItem* editted = new QStandardItem(v.toString());
            m_treeNodeVars->insertRow(item.row() + 1, editted);
            m_treeNodeVars->removeRow(item.row());
        }
    }
    else if (parent.row() == m_csRow)
    {

    }
    else if (parent.row() == m_matRow)
    {
//        MaterialDialog* dialog = new MaterialDialog(this->materials.at(item.row()),
//                                                    this->matDict, this->varDict, this);
//        if (dialog->exec() == QDialog::Accepted)
//        {
//            Material m = dialog->data();
//            this->materials[item.row()] = m;
//            this->matDict[m.name()] = m;
//            QStandardItem* editted = new QStandardItem(m.toString());
//            treeNodeMats->insertRow(item.row() + 1, editted);
//            treeNodeMats->removeRow(item.row());
//        }
    }
    else {
        qWarning("%s", QString(tr("Unknown item in declarations!")).toStdString().c_str());
    }
}

void DeclarationsWidget::on_DeclarationTree_doubleClicked(const QModelIndex &index)
{
    this->createEditDialog(index);
}
