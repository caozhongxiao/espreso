#include "regionpairdialog.h"
#include "ui_regionpairdialog.h"

#include <QLabel>
#include <QComboBox>
#include <QFormLayout>

#include "../declarations/datatypeeditwidget.h"

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

RegionPairDialog::~RegionPairDialog()
{
    delete ui;
}

RegionPairDialog* RegionPairDialog::createRegionMaterial(ECFObject* map, Mesh* mesh, ECFObject* materials)
{
    return new RegionPairDialog(ECFDataType::MATERIAL, map, mesh, materials);
}

RegionPairDialog* RegionPairDialog::createRegionExpression(ECFObject* map, Mesh* mesh)
{
    return new RegionPairDialog(ECFDataType::EXPRESSION, map, mesh, nullptr);
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
        value->setValue(w->value().toStdString());
    }

    this->m_region = QString::fromStdString(key);

    QDialog::accept();
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
