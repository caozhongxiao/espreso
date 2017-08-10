#include "materialdialog.h"
#include "ui_materialdialog.h"

#include "../../data/common.h"
#include "datatypewidget.h"
#include <QStringList>

MaterialDialog::MaterialDialog(const QHash<QString, Material>& matDict, const QHash<QString, Variable>& varDict, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MaterialDialog)
{
    ui->setupUi(this);

    this->matDict = matDict;
    this->varDict = varDict;

    QRegExpValidator* nameValidator = new QRegExpValidator(QRegExp(REGEXPR_NAME), this);
    ui->editMatName->setValidator(nameValidator);

    this->setupProperties();

    ui->frProperty->layout()->addWidget(new DataTypeWidget(this->varDict, ui->frProperty));

}

MaterialDialog::MaterialDialog(const Material& material, const QHash<QString, Material>& matDict, const QHash<QString, Variable>& varDict, QWidget *parent) :
    MaterialDialog(matDict, varDict, parent)
{
    this->ui->editMatName->setText(material.name());
    this->ui->editMatDesc->setPlainText(material.description());
}

MaterialDialog::~MaterialDialog()
{
    delete ui;
}

Material MaterialDialog::data() const
{
    Material m(ui->editMatName->text(), ui->editMatDesc->toPlainText());

    return m;
}

void MaterialDialog::on_btnPropAdd_pressed()
{
    QList<QStandardItem*> list;
    QStandardItem* cell1 = new QStandardItem("name");
    QStandardItem* cell2 = new QStandardItem("unit");
    QStandardItem* cell3 = new QStandardItem("model");
    QStandardItem* cell4 = new QStandardItem("button");
    list << cell1 << cell2 << cell3 << cell4;
    this->tableModel->appendRow(list);
}

void MaterialDialog::on_btnPropDel_pressed()
{
    QModelIndexList indices = ui->tableProperties->selectionModel()->selectedIndexes();

    foreach (QModelIndex index, indices) {
        this->tableModel->removeRow(index.row());
        this->deletedProperties.insert(index.row(), index.row());
    }
}

void MaterialDialog::setupProperties()
{
    this->tableModel = new QStandardItemModel(this);
    ui->tableProperties->setModel(this->tableModel);
    QStringList tableHeader;
    tableHeader << "Name" << "Unit" << "Model" << "Edit";
    this->tableModel->setHorizontalHeaderLabels(tableHeader);

    ui->tableProperties->setEditTriggers(QAbstractItemView::NoEditTriggers);

    IsotropicProperty* thermal = new IsotropicProperty("Thermal Conductivity", "W*m^-1*K^-1", new IsotropicMatrix());
    BasicProperty* density = new BasicProperty("Density", "kg*m^-3", new BasicMatrix());
    BasicProperty* cp = new BasicProperty("Specific Heat", "J*kg^-1*K^-1", new BasicMatrix());

    this->properties.append(thermal);
    this->properties.append(density);
    this->properties.append(cp);

    foreach (MaterialProperty* mp, this->properties) {
        QList<QStandardItem*> list;
        QStandardItem* cell1 = new QStandardItem(mp->name());
        QStandardItem* cell2 = new QStandardItem(mp->unit());
        QStandardItem* cell3 = new QStandardItem("Model");
        list << cell1 << cell2 << cell3;
        this->tableModel->appendRow(list);
        QPushButton* btn = new QPushButton("Edit", ui->tableProperties);
        int row = this->tableModel->rowCount() - 1;
        connect(btn, &QPushButton::pressed,
                [=] () { tableBtnPressed(row); });
        QModelIndex btnIndex = this->tableModel->index(row, 3);
        ui->tableProperties->setIndexWidget(btnIndex, btn);
    }
}

void MaterialDialog::tableBtnPressed(int row)
{
    MaterialProperty* p = this->properties.at(row);
    p->accept(this);
}

void MaterialDialog::visit(const BasicProperty& property)
{

}

void MaterialDialog::visit(const IsotropicProperty& property)
{

}

void MaterialDialog::visit(const DiagonalProperty& property)
{

}

void MaterialDialog::visit(const SymmetricProperty& property)
{

}

void MaterialDialog::visit(const AnisotropicProperty& property)
{

}
