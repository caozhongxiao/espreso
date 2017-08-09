#include "datatypewidget.h"
#include "ui_datatypewidget.h"
#include "doubletabledelegate.h"

#include <QMessageBox>
#include "../../expression.h"


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

        if (!Expression::isValid(ui->editFunction->text().toStdString(), Common::fnVariables()))
        {
            QMessageBox::warning(this, tr("Error"), tr("Incorrect format of function formula!"));
            return -1;
        }
    }
    else if (ui->cmbType->currentIndex() == 3)
    {
        int rowCount = this->piecewiseModel->rowCount();
        for (int row = 0; row < rowCount; ++row)
        {
            QString fn = this->piecewiseModel->item(row, 2)->text();
            if (!Expression::isValid(fn.toStdString(), Common::fnVariables()))
            {
                QMessageBox::warning(this, tr("Error"), tr("Incorrect format of function formula in table!"));
                return -1;
            }
        }
    }

    return 0;
}

void DataTypeWidget::setData(const DataType* data)
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

void DataTypeWidget::setupTableData(const DataType* data)
{
    const TableType* table = (const TableType*)data;
    for (auto it = table->data().cbegin(); it != table->data().cend(); ++it) {
        QList<QStandardItem*> list;
        QStandardItem* cell1 = new QStandardItem(it->first);
        QStandardItem* cell2 = new QStandardItem(it->second);
        list << cell1 << cell2;
        this->tableModel->appendRow(list);
    }
}

DataType* DataTypeWidget::collectTableData() const
{
    int rowCount = this->tableModel->rowCount();
    QVector<QPair<QString, QString> > rowVector;
    for (int row = 0; row < rowCount; ++row)
    {
        QPair<QString, QString> tuple;
        tuple.first = this->tableModel->item(row, 0)->text();
        tuple.second = this->tableModel->item(row, 1)->text();
        rowVector.append(tuple);
    }

    return new TableType(rowVector);
}

void DataTypeWidget::setupPiecewiseData(const DataType* data)
{
    const PiecewiseFunctionType* table = (const PiecewiseFunctionType*)data;
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
    QVector<QVector<QString> > rowVector;
    for (int row = 0; row < rowCount; ++row)
    {
        QVector<QString> triple;
        triple.append(this->piecewiseModel->item(row, 0)->text());
        triple.append(this->piecewiseModel->item(row, 1)->text());
        triple.append(this->piecewiseModel->item(row, 2)->text());
        rowVector.append(triple);
    }

    return new PiecewiseFunctionType(rowVector);
}

void DataTypeWidget::setupVariableLinkData(const DataType* data)
{
    const VariableLinkType* varLink = (const VariableLinkType*)data;
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

void DataTypeWidget::visit(const ConstantType& type)
{
    this->ui->cmbType->setCurrentIndex(0);
}

void DataTypeWidget::visit(const FunctionType& type)
{
    this->ui->cmbType->setCurrentIndex(1);
}

void DataTypeWidget::visit(const TableType& type)
{
    this->ui->cmbType->setCurrentIndex(2);
}

void DataTypeWidget::visit(const PiecewiseFunctionType& type)
{
    this->ui->cmbType->setCurrentIndex(3);
}

void DataTypeWidget::visit(const VariableLinkType& type)
{
    this->ui->cmbType->setCurrentIndex(4);
}
