#include "datatypewidget.h"
#include "ui_datatypewidget.h"
#include "doubletabledelegate.h"

#include <QMessageBox>


DataTypeWidget::DataTypeWidget(const QHash<QString, Variable>& varDict, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DataTypeWidget)
{
    ui->setupUi(this);

    this->varDict = varDict;

    this->setupCheckbox();

    // SETUP GUI ELEMENTS OF INDIVIDUAL DATA TYPES (COMBOBOX CHOICES)
    QRegExpValidator* constantValidator = new QRegExpValidator(QRegExp(REGEXPR_DOUBLE), this);
    ui->editConstant->setValidator(constantValidator);

    this->tableModel = new QStandardItemModel(this);
    ui->tblTable->setModel(this->tableModel);
    ui->tblTable->setItemDelegate(new DoubleTableDelegate());
    QStringList tableHeader;
    tableHeader << "x" << "f(x)";
    this->tableModel->setHorizontalHeaderLabels(tableHeader);

    this->piecewiseModel = new QStandardItemModel(this);
    ui->tblPiecewise->setModel(this->piecewiseModel);
    QStringList piecewiseHeader;
    piecewiseHeader << tr("lower bound") << tr("upper bound") << tr("function");
    this->piecewiseModel->setHorizontalHeaderLabels(piecewiseHeader);
    ui->tblPiecewise->setItemDelegateForColumn(0, new DoubleTableDelegate());
    ui->tblPiecewise->setItemDelegateForColumn(1, new DoubleTableDelegate());

    foreach (QString key, this->varDict.keys()) {
        ui->cmbVariable->addItem(key);
    }
    if (!this->varDict.isEmpty())
        ui->cmbVariable->setCurrentIndex(0);
}

DataTypeWidget::DataTypeWidget(DataType* data, const QHash<QString, Variable>& varDict, QWidget *parent) :
    DataTypeWidget(varDict, parent)
{
    data->accept(this);
    this->setData(data);
}

DataTypeWidget::~DataTypeWidget()
{
    delete ui;
}

void DataTypeWidget::setupCheckbox()
{
    ui->cmbType->addItem(tr("Constant"));
    ui->cmbType->addItem(tr("Function"));
    ui->cmbType->addItem(tr("Table"));
    ui->cmbType->addItem(tr("Piecewise function"));
    if (!this->varDict.isEmpty())
        ui->cmbType->addItem(tr("Variable"));
    ui->cmbType->setCurrentIndex(0);
}

void DataTypeWidget::setDataType(DataType* dataType)
{
    dataType->accept(this);
    this->setData(dataType);
}

DataType* DataTypeWidget::data() const
{
    DataType* data;
    switch(ui->cmbType->currentIndex())
    {
        case 0:
            data = new ConstantType(ui->editConstant->text());
            break;
        case 1:
            data = new FunctionType(ui->editFunction->text());
            break;
        case 2:
            data = this->collectTableData();
            break;
        case 3:
            data = this->collectPiecewiseData();
            break;
        case 4:
            data = new VariableLinkType(ui->cmbVariable->itemText(ui->cmbVariable->currentIndex()));
            break;
        default:
            qWarning("%s", QString(tr("Dialog: Cannot retrieve data from an unknown type!")).toStdString().c_str());
    }

    return data;
}

DataType* DataTypeWidget::checkedData()
{
    if (this->check() == -1)
        return nullptr;
    return this->data();
}

void DataTypeWidget::on_cmbType_currentIndexChanged(int index)
{
    switch(index)
    {
        case 0:
            ui->editFunction->hide();
            ui->frameTable->hide();
            ui->framePiecewise->hide();
            ui->cmbVariable->hide();
            ui->editConstant->clear();
            ui->editConstant->show();
            break;
        case 1:
            ui->editConstant->hide();
            ui->frameTable->hide();
            ui->framePiecewise->hide();
            ui->cmbVariable->hide();
            ui->editFunction->clear();
            ui->editFunction->show();
            break;
        case 2:
            ui->editFunction->hide();
            ui->editConstant->hide();
            ui->framePiecewise->hide();
            ui->cmbVariable->hide();
            this->tableModel->removeRows(0, tableModel->rowCount());
            ui->frameTable->show();
            break;
        case 3:
            ui->editFunction->hide();
            ui->editConstant->hide();
            ui->frameTable->hide();
            ui->cmbVariable->hide();
            this->piecewiseModel->removeRows(0, piecewiseModel->rowCount());
            ui->framePiecewise->show();
            break;
        case 4:
            ui->editFunction->hide();
            ui->editConstant->hide();
            ui->frameTable->hide();
            ui->framePiecewise->hide();
            ui->cmbVariable->show();
            break;
        default:
            qWarning("%s (%d)", QString(tr("Dialog: Unknown type selected!")).toStdString().c_str(), index);
    }
}

int DataTypeWidget::check()
{
    // FINAL CHECK BEFORE FORM SUBMISSION
    if (ui->cmbType->currentIndex() == 1)
    {
        if (ui->editFunction->text().isEmpty())
        {
            QMessageBox::warning(this, tr("Error"), tr("Function formula not entered!"));
            return -1;
        }

        //TODO Expression validation
    }
    else if (ui->cmbType->currentIndex() == 3)
    {
        int rowCount = this->piecewiseModel->rowCount();
        for (int row = 0; row < rowCount; ++row)
        {
            //TODO Expression validation
        }
    }

    return 1;
}

void DataTypeWidget::setData(DataType* data)
{
    // FILL THE DIALOG WITH DATA OF CURRENT ITEM
    switch(ui->cmbType->currentIndex())
    {
        case 0:
            ui->editConstant->setText(data->toString());
            break;
        case 1:
            ui->editFunction->setText(data->toString());
            break;
        case 2:
            this->setupTableData(data);
            break;
        case 3:
            this->setupPiecewiseData(data);
            break;
        case 4:
            this->setupVariableLinkData(data);
            break;
        default:
            qWarning("%s", QString(tr("Dialog: Cannot display data of an unknown type!")).toStdString().c_str());
    }
}

void DataTypeWidget::setupTableData(DataType* data)
{
    this->tableModel->clear();
    TableType* table = (TableType*)data;
    for (auto it = table->data().cbegin(); it != table->data().cend(); ++it) {
        QList<QStandardItem*> list;
        QStandardItem* cell1 = new QStandardItem(it->at(0));
        QStandardItem* cell2 = new QStandardItem(it->at(1));
        list << cell1 << cell2;
        this->tableModel->appendRow(list);
    }
}

DataType* DataTypeWidget::collectTableData() const
{
    int rowCount = this->tableModel->rowCount();
    QList<QList<QString> > rowVector;
    for (int row = 0; row < rowCount; ++row)
    {
        QList<QString> list;
        list << this->tableModel->item(row, 0)->text();
        list << this->tableModel->item(row, 1)->text();
        rowVector.append(list);
    }

    return new TableType(rowVector);
}

void DataTypeWidget::setupPiecewiseData(DataType* data)
{
    this->piecewiseModel->clear();
    PiecewiseFunctionType* table = (PiecewiseFunctionType*)data;
    for (auto it = table->data().cbegin(); it != table->data().cend(); ++it) {
        QList<QStandardItem*> list;
        QStandardItem* cell1 = new QStandardItem(it->at(0));
        QStandardItem* cell2 = new QStandardItem(it->at(1));
        QStandardItem* cell3 = new QStandardItem(it->at(2));
        list << cell1 << cell2 << cell3;
        this->piecewiseModel->appendRow(list);
    }
}

DataType* DataTypeWidget::collectPiecewiseData() const
{
    int rowCount = this->piecewiseModel->rowCount();
    QList<QList<QString> > rowVector;
    for (int row = 0; row < rowCount; ++row)
    {
        QList<QString> triple;
        triple.append(this->piecewiseModel->item(row, 0)->text());
        triple.append(this->piecewiseModel->item(row, 1)->text());
        triple.append(this->piecewiseModel->item(row, 2)->text());
        rowVector.append(triple);
    }

    return new PiecewiseFunctionType(rowVector);
}

void DataTypeWidget::setupVariableLinkData(DataType* data)
{
    VariableLinkType* varLink = (VariableLinkType*)data;
    QList<QString> keys = this->varDict.keys();
    ui->cmbVariable->setCurrentIndex(keys.indexOf(varLink->toString()));
}

void DataTypeWidget::on_btnTableAdd_pressed()
{
    QList<QStandardItem*> list;
    QStandardItem* cell1 = new QStandardItem("0");
    QStandardItem* cell2 = new QStandardItem("0");
    list << cell1 << cell2;
    this->tableModel->appendRow(list);
}

void DataTypeWidget::on_btnTableDel_pressed()
{
    QModelIndexList indices = ui->tblTable->selectionModel()->selectedIndexes();

    foreach (QModelIndex index, indices) {
        this->tableModel->removeRow(index.row());
    }
}

void DataTypeWidget::on_btnPiecewiseAdd_pressed()
{
    QList<QStandardItem*> list;
    QStandardItem* cell1 = new QStandardItem("0");
    QStandardItem* cell2 = new QStandardItem("1");
    QStandardItem* cell3 = new QStandardItem("x");
    list << cell1 << cell2 << cell3;
    this->piecewiseModel->appendRow(list);
}

void DataTypeWidget::on_btnPiecewiseDel_pressed()
{
    QModelIndexList indices = ui->tblPiecewise->selectionModel()->selectedIndexes();

    foreach (QModelIndex index, indices) {
        this->piecewiseModel->removeRow(index.row());
    }
}

void DataTypeWidget::visit(ConstantType& type)
{
    this->ui->cmbType->setCurrentIndex(0);
}

void DataTypeWidget::visit(FunctionType& type)
{
    this->ui->cmbType->setCurrentIndex(1);
}

void DataTypeWidget::visit(TableType& type)
{
    this->ui->cmbType->setCurrentIndex(2);
}

void DataTypeWidget::visit(PiecewiseFunctionType& type)
{
    this->ui->cmbType->setCurrentIndex(3);
}

void DataTypeWidget::visit(VariableLinkType& type)
{
    this->ui->cmbType->setCurrentIndex(4);
}
