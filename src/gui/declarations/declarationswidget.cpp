#include "declarationswidget.h"
#include "ui_declarationswidget.h"

#include "variabledialog.h"
#include "../data/variable.h"
#include "material/materialdialog.h"
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

void DeclarationsWidget::initialize(PhysicsConfiguration *physics)
{
    this->m_physics = physics;
    this->m_materials = &m_physics->materials;
    for (auto it = this->m_materials->begin(); it != this->m_materials->end(); ++it)
    {
        bool ok;
        int mid = QString::fromStdString(it->first).toInt(&ok);
        if (ok && mid > this->m_materialID) this->m_materialID = mid;
        this->m_materialNames.append(it->second.getParameter("name")->getValue());
        this->m_materialIDs.append(it->first);
        QStandardItem* item = new QStandardItem(
                    QString::fromStdString(
                        it->second.getParameter("name")->getValue()
                        )
                    );
        m_treeNodeMats->appendRow(item);
    }
    this->m_materialID++;

    this->m_initialized = true;
}

void DeclarationsWidget::setPhysics(PhysicsConfiguration *physics)
{
    if (!this->m_initialized)
    {
        qFatal("DeclarationsWidget: Method initialize should be called before one can change physics.");
        return;
    }

    this->m_physics = physics;
    physics->materials.clear();

    for (auto it = m_materials->begin(); it != m_materials->end(); ++it)
    {
        physics->materials[it->first] = it->second;
    }

    this->m_materials = &m_physics->materials;
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
        MaterialConfiguration* material = this->createMaterial();
        MaterialDialog* dialog = new MaterialDialog(material, m_materialNames, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            QStandardItem* item = new QStandardItem(
                        QString::fromStdString(
                            material->getParameter(std::string("name"))->getValue()
                            )
                        );
            m_treeNodeMats->appendRow(item);
            this->m_materialNames.append(
                        this->toUpper(
                            material->getParameter(std::string("name"))->getValue()
                            )
                        );
        }
        else
        {
            this->removeMaterial(m_materialIDs.size() - 1);
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
        this->removeMaterial(clicked.row());
        m_treeNodeMats->removeRow(clicked.row());
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
        std::string mID = m_materialIDs.at(item.row());
        m_materialNames.remove(item.row());
        MaterialConfiguration* material = &(*m_materials)[mID];
        MaterialDialog* dialog = new MaterialDialog(material, m_materialNames, this);
        if (dialog->exec() == QDialog::Accepted)
        {
            QStandardItem* editted = new QStandardItem(
                        QString::fromStdString(
                            material->getParameter(std::string("name"))->getValue()
                            )
                        );
            m_treeNodeMats->insertRow(item.row() + 1, editted);
            m_treeNodeMats->removeRow(item.row());
        }

        m_materialNames.insert(
                    item.row(),
                    this->toUpper(
                        material->getParameter(std::string("name"))->getValue()
                        ));
    }
    else {
        qWarning("%s", QString(tr("Unknown item in declarations!")).toStdString().c_str());
    }
}

void DeclarationsWidget::on_DeclarationTree_doubleClicked(const QModelIndex &index)
{
    this->createEditDialog(index);
}

MaterialConfiguration* DeclarationsWidget::createMaterial()
{
    std::string mID = this->createMaterialId();

    (*this->m_materials)[mID];// = MaterialConfiguration();

//    if (this->m_physics == nullptr)
//    {
//        (*this->m_materials)[mID];// = MaterialConfiguration();
//    }
//    else
//    {
//        qCritical("NOT SUPPORTED YET!");
////        (*this->m_materials)[mID] = MaterialConfiguration(
////                    this->m_physics->dimension,
////                    this->m_physics->physical_model
////                    );
//    }
    this->m_materialIDs.append(mID);

    return &(*this->m_materials)[mID];
}

std::string DeclarationsWidget::createMaterialId()
{
    return std::to_string(m_materialID++);
}

void DeclarationsWidget::removeMaterial(int index)
{
    std::string key = this->m_materialIDs.at(index);
    this->m_materials->erase(key);
    this->m_materialIDs.remove(index);
    this->m_materialNames.remove(index);
}

std::string DeclarationsWidget::toUpper(const std::string& text)
{
    return QString::fromStdString(text)
            .toUpper()
            .toStdString();
}
