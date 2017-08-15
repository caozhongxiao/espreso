#include "materialdialog.h"
#include "ui_materialdialog.h"

#include "../../data/common.h"
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

    this->basicPropertyWidget = new DataTypeWidget(this->varDict, ui->frPropertyDetail);
    this->matrixPropertyWidget = new MatrixPropertyWidget(this->varDict, ui->frPropertyDetail);
    ui->frPropertyDetail->layout()->addWidget(this->basicPropertyWidget);
    ui->frPropertyDetail->layout()->addWidget(this->matrixPropertyWidget);
    ui->frProperty->hide();
    this->basicPropertyWidget->hide();
    this->matrixPropertyWidget->hide();
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

    IsotropicModel im;
    im.kxx = new ConstantType("0");
    BasicModel bm;
    bm.kxx = new ConstantType("0");
    DiagonalModel dm;
    dm.kxx = new ConstantType("0");
    dm.kyy = new ConstantType("0");
    dm.kzz = new ConstantType("0");
    MaterialProperty* thermal = new IsotropicProperty("Thermal Conductivity", "W*m^-1*K^-1", im);
    MaterialProperty* thermal2 = new DiagonalProperty("Thermal Conductivity 2", "W*m^-1*K^-1", dm);
    MaterialProperty* density = new BasicProperty("Density", "kg*m^-3", bm);
    MaterialProperty* cp = new BasicProperty("Specific Heat", "J*kg^-1*K^-1", bm);
    this->properties.append(thermal);
    this->properties.append(thermal2);
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
    this->activeProperty = row;
    this->activePropertyWidget = 0;
    p->accept(this);
}

void MaterialDialog::visit(BasicProperty& p)
{
    this->basicPropertyWidget->setDataType(p.model().kxx);
    this->matrixPropertyWidget->hide();
    ui->frProperty->show();
    this->basicPropertyWidget->show();
}

void MaterialDialog::visit(IsotropicProperty& p)
{
    this->matrixPropertyDetail(&p);
}

void MaterialDialog::visit(DiagonalProperty& p)
{
    this->matrixPropertyDetail(&p);
}

void MaterialDialog::visit(SymmetricProperty& p)
{
    this->matrixPropertyDetail(&p);
}

void MaterialDialog::visit(AnisotropicProperty& p)
{
    this->matrixPropertyDetail(&p);
}

void MaterialDialog::matrixPropertyDetail(MaterialProperty *p)
{
    this->activePropertyWidget = 1;
    this->basicPropertyWidget->hide();
    this->matrixPropertyWidget->setData(p);
    this->matrixPropertyWidget->show();
    ui->frProperty->show();
}

void MaterialDialog::on_btnPropSave_pressed()
{
    int ret = 1;
    switch(activePropertyWidget)
    {
        case 0:
            ret = this->serveBasicPropertyWidget();
            break;
        case 1:
            ret = this->serveMatrixPropertyWidget();
            break;
        default:
            qWarning("%s", QString(tr("Property detail: Unknown property widget!")).toStdString().c_str());
    }
    if (ret)
        ui->frProperty->hide();
}

int MaterialDialog::serveBasicPropertyWidget()
{
    if (this->basicPropertyWidget->check() == -1)
        return 0;
    DataType* result = this->basicPropertyWidget->data();
    MaterialProperty* p = this->properties.at(this->activeProperty);
    QVector<DataType*> data;
    data.append(result);
    p->setModelData(data);
    return 1;
}

int MaterialDialog::serveMatrixPropertyWidget()
{
    MaterialProperty* result = this->matrixPropertyWidget->data();
    this->properties[this->activeProperty] = result;
    qInfo("%s", result->name().toStdString().c_str());
    return 1;
}
