#include "regionpairdialog.h"
#include "ui_regionpairdialog.h"

#include <QLabel>
#include <QComboBox>
#include <QFormLayout>
#include <QMessageBox>
#include <QDebug>
#include <QPushButton>

#include "../declarations/datatypeeditwidget.h"
#include "regionobjectwidget.h"

using namespace espreso;

RegionPairDialog::RegionPairDialog(ECFDataType value, ECFObject* map,
                                 Mesh* mesh, ECFObject* scope, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RegionPairDialog)
{
    ui->setupUi(this);

    this->m_first = ECFDataType::REGION;
    this->m_second = value;

    this->m_map = map;
    this->m_scope = scope;
    this->m_mesh = mesh;

    this->m_first_widget = this->uiValue(m_first, ui->first);
    this->m_second_widget = this->uiValue(m_second, ui->second);
}

RegionPairDialog::RegionPairDialog(ECFParameter* pair, ECFDataType value,
                         ECFObject* map, Mesh* mesh,
                         ECFObject* scope, QWidget *parent)
    : RegionPairDialog(value, map, mesh, scope, parent)
{
    this->m_first_widget->setEnabled(false);
    QComboBox* cmb = static_cast<QComboBox*>(this->m_first_widget);
    cmb->clear();
    cmb->addItem(QString::fromStdString(pair->name));

    if (this->m_second == ECFDataType::MATERIAL)
    {
        QComboBox* mat = static_cast<QComboBox*>(this->m_second_widget);
        std::string value = pair->getValue();
        int i = 0;
        for (auto it = scope->parameters.begin(); it != scope->parameters.end(); ++it)
        {
            if ((*it)->getParameter("name")->getValue().compare(value) == 0)
            {
                mat->setCurrentIndex(i);
                break;
            }
            i++;
        }
    }
    else if (this->m_second == ECFDataType::EXPRESSION)
    {
        DataTypeEditWidget* expr = static_cast<DataTypeEditWidget*>(this->m_second_widget);
        std::string value = pair->getValue();
        expr->setValue(QString::fromStdString(value));
    }
}

RegionPairDialog::RegionPairDialog(ECFObject *map, Mesh *mesh, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RegionPairDialog)
{
    ui->setupUi(this);

    this->m_first = ECFDataType::REGION;
    this->m_second = ECFDataType::LOAD_STEP;

    this->m_map = map;
    this->m_mesh = mesh;

    this->m_first_widget = this->uiValue(m_first, ui->first);

    ECFObject* obj = static_cast<ECFObject*>(this->m_map->getPattern());
    RegionObjectWidget* row = new RegionObjectWidget(obj);
    row->init();
    this->m_object = obj;
    this->m_second_widget = row;
    ui->second->addWidget(row);
}

RegionPairDialog::RegionPairDialog(ECFObject *pair, ECFObject *map, Mesh *mesh, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RegionPairDialog)
{
    ui->setupUi(this);

    this->m_first = ECFDataType::REGION;
    this->m_second = ECFDataType::LOAD_STEP;

    this->m_map = map;
    this->m_mesh = mesh;

    this->m_first_widget = this->uiValue(m_first, ui->first);
    this->m_first_widget->setEnabled(false);

    RegionObjectWidget* row = new RegionObjectWidget(pair);
    row->init();
    this->m_second_widget = row;
    ui->second->addWidget(row);
}


RegionPairDialog::~RegionPairDialog()
{
    delete ui;
}

RegionPairDialog* RegionPairDialog::createRegionMaterial(ECFObject* map, Mesh* mesh, ECFObject* materials, ECFParameter* pair)
{
    if (pair == nullptr)
        return new RegionPairDialog(ECFDataType::MATERIAL, map, mesh, materials);
    else
        return new RegionPairDialog(pair, ECFDataType::MATERIAL, map, mesh, materials);
}

RegionPairDialog* RegionPairDialog::createRegionExpression(ECFObject* map, Mesh* mesh, ECFParameter* pair)
{
    if (pair == nullptr)
        return new RegionPairDialog(ECFDataType::EXPRESSION, map, mesh, nullptr);
    else
        return new RegionPairDialog(pair, ECFDataType::EXPRESSION, map, mesh, nullptr);
}

RegionPairDialog* RegionPairDialog::createRegionObject(ECFObject *map, Mesh *mesh, ECFObject *pair)
{
    if (pair == nullptr)
        return new RegionPairDialog(map, mesh);
    else
        return new RegionPairDialog(pair, map, mesh);
}

QWidget* RegionPairDialog::uiValue(ECFDataType type, QLayout* layout)
{
    QWidget* ret;

    if (type == ECFDataType::REGION)
    {
        layout->addWidget(new QLabel(tr("Region:"), this));
        QComboBox* cmb = new QComboBox(this);
        for (auto it = m_mesh->regions().begin(); it != m_mesh->regions().end(); it++)
        {
            bool found = false;
            for (auto key = this->m_map->parameters.begin(); key != this->m_map->parameters.end(); key++)
            {
                if ((*key)->name.compare((*it)->name) == 0)
                {
                    found = true;
                    break;
                }
            }
            if (found) continue;
            cmb->addItem(QString::fromStdString((*it)->name));
        }
        layout->addWidget(cmb);
        ret = cmb;
    }

    if (type == ECFDataType::MATERIAL)
    {
        layout->addWidget(new QLabel(tr("Material:"), this));
        QComboBox* cmb = new QComboBox(this);
        for (auto it = this->m_scope->parameters.begin(); it != this->m_scope->parameters.end(); ++it)
        {
            cmb->addItem(QString::fromStdString((*it)->getParameter("name")->getValue()));
        }
        if (this->m_scope->parameters.size() > 0) cmb->setCurrentIndex(0);
        else ui->buttonBox->setEnabled(false);
        layout->addWidget(cmb);
        ret = cmb;
    }

    if (type == ECFDataType::EXPRESSION)
    {
        ECFObject* tmp = static_cast<ECFObject*>(this->m_map->getParameter("***"));
        DataTypeEditWidget* w = new DataTypeEditWidget(tmp->metadata.variables, this);
        this->m_map->dropParameter(tmp);
        QWidget* container = new QWidget;
        QFormLayout* fl = new QFormLayout;
        fl->addRow(tr("Type:"), w->createComboBox(this));
        fl->addRow(tr("Value:"), w);
        container->setLayout(fl);
        layout->addWidget(container);
        ret = w;
    }

    return ret;
}

void RegionPairDialog::accept()
{
    std::string key = this->getKey();

    ECFParameter* value = this->m_map->getParameter(key);

    if (this->m_second == ECFDataType::MATERIAL)
    {
        QComboBox* cmb = static_cast<QComboBox*>(this->m_second_widget);
        int index = cmb->currentIndex();
        value->setValue(this->m_scope->parameters[index]->name);
    }

    if (this->m_second == ECFDataType::EXPRESSION)
    {
        DataTypeEditWidget* w = static_cast<DataTypeEditWidget*>(this->m_second_widget);
        if (!w->isValid())
        {
            if (value->getValue().empty()) this->m_map->dropParameter(value);
            this->displayError(w->errorMessage());
            return;
        }
        value->setValue(w->value().toStdString());
    }

    if (this->m_second == ECFDataType::LOAD_STEP)
    {
        RegionObjectWidget* row = static_cast<RegionObjectWidget*>(this->m_second_widget);
        if (!row->isValid())
        {
            if (this->m_object != nullptr) this->m_map->dropParameter(value);
            this->displayError(row->errorMessage());
            return;
        }
        row->save();
        if (this->m_object != nullptr)
        {
            ECFObject* newObj = static_cast<ECFObject*>(value);
            int pi = 0;
            for (auto p = this->m_object->parameters.begin(); p != this->m_object->parameters.end(); ++p)
            {
                newObj->parameters[pi++] = (*p);
            }
        }
    }

    this->m_region = QString::fromStdString(key);

    QDialog::accept();
}

void RegionPairDialog::displayError(const QString& msg)
{
    QMessageBox msgbox;
    msgbox.setWindowTitle(tr("Error"));
    msgbox.setText(msg);
    msgbox.exec();
}

QString RegionPairDialog::region()
{
    return this->m_region;
}

std::string RegionPairDialog::getKey()
{
    if (this->m_first == ECFDataType::REGION)
    {
        QComboBox* cmb = static_cast<QComboBox*>(this->m_first_widget);
        return cmb->currentText().toStdString();
    }

    return "";
}
