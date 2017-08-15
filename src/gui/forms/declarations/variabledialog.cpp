#include "variabledialog.h"
#include "ui_variabledialog.h"

#include "../../data/common.h"
#include <QMessageBox>
#include "../../expression.h"
#include "doubletabledelegate.h"
#include <QPair>

VariableDialog::VariableDialog(const QHash<QString, Variable>& varDict, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VariableDialog)
{
    ui->setupUi(this);

    // SETUP ITEMS
    ui->cmbType->addItem(tr("Constant"));
    ui->cmbType->addItem(tr("Function"));
    ui->cmbType->addItem(tr("Table"));
    ui->cmbType->addItem(tr("Piecewise function"));
    ui->cmbType->setCurrentIndex(0);

    // SETUP GUI ELEMENTS OF INDIVIDUAL DATA TYPES (COMBOBOX CHOICES)
    QRegExpValidator* varNameValidator = new QRegExpValidator(QRegExp(REGEXPR_VARNAME), this);
    ui->editName->setValidator(varNameValidator);

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

    this->varDict = varDict;

    fnVars.push_back("x");
    fnVars.push_back("y");
    fnVars.push_back("z");
    fnVars.push_back("time");
    fnVars.push_back("temperature");
}

VariableDialog::VariableDialog(const Variable& var, const QHash<QString, Variable>& varDict, QWidget *parent) :
    VariableDialog(varDict, parent)
{
    this->varDict.remove(var.name());
    ui->editName->setText(var.name());
    var.data()->accept(this);
    this->setData(var);
}

//VariableDialog::VariableDialog(QVariant data, QWidget *parent) :
//    VariableDialog(parent)
//{
//    if (!data.canConvert<Variable>())
//    {
//        qWarning("%s", tr("VariableDialog: No data provided for the dialog!").toStdString().c_str());
//        return;
//    }
//    Variable var = data.value<Variable>();
//    ui->editName->setText(var.name());
//    ui->editConstant->setText(var.data()->toString());
//}

VariableDialog::~VariableDialog()
{
    delete ui;
}

Variable VariableDialog::data()
{
    QString name = ui->editName->text();
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
        default:
            qWarning("%s", QString(tr("Dialog: Cannot retrieve data from an unknown type!")).toStdString().c_str());
    }

    return Variable(name, data);
}

void VariableDialog::on_cmbType_currentIndexChanged(int index)
{
    switch(index)
    {
        case 0:
            ui->editFunction->hide();
            ui->frameTable->hide();
            ui->framePiecewise->hide();
            ui->editConstant->clear();
            ui->editConstant->show();
            break;
        case 1:
            ui->editConstant->hide();
            ui->frameTable->hide();
            ui->framePiecewise->hide();
            ui->editFunction->clear();
            ui->editFunction->show();
            break;
        case 2:
            ui->editFunction->hide();
            ui->editConstant->hide();
            ui->framePiecewise->hide();
            this->tableModel->removeRows(0, tableModel->rowCount());
            ui->frameTable->show();
            break;
        case 3:
            ui->editFunction->hide();
            ui->editConstant->hide();
            ui->frameTable->hide();
            this->piecewiseModel->removeRows(0, piecewiseModel->rowCount());
            ui->framePiecewise->show();
            break;
        default:
            qWarning("%s (%d)", QString(tr("Dialog: Unknown type selected!")).toStdString().c_str(), index);
    }
}

void VariableDialog::accept()
{
    // FINAL CHECK BEFORE FORM SUBMISSION
    if (this->varDict.contains(ui->editName->text()))
    {
        QMessageBox::warning(this, tr("Error"), tr("Chosen variable name is already in use!"));
        return;
    }
    else if (ui->cmbType->currentIndex() == 1)
    {
        if (ui->editFunction->text().isEmpty())
        {
            QMessageBox::warning(this, tr("Error"), tr("Function formula not entered!"));
            return;
        }

        if (!Expression::isValid(ui->editFunction->text().toStdString(), this->fnVars))
        {
            QMessageBox::warning(this, tr("Error"), tr("Incorrect format of function formula!"));
            return;
        }
    }
    else if (ui->cmbType->currentIndex() == 3)
    {
        int rowCount = this->piecewiseModel->rowCount();
        for (int row = 0; row < rowCount; ++row)
        {
            QString fn = this->piecewiseModel->item(row, 2)->text();
            if (!Expression::isValid(fn.toStdString(), this->fnVars))
            {
                QMessageBox::warning(this, tr("Error"), tr("Incorrect format of function formula in table!"));
                return;
            }
        }
    }

    QDialog::accept();
}

void VariableDialog::setData(const Variable &var)
{
    // FILL THE DIALOG WITH DATA OF CURRENT ITEM
    switch(ui->cmbType->currentIndex())
    {
        case 0:
            ui->editConstant->setText(var.data()->toString());
            break;
        case 1:
            ui->editFunction->setText(var.data()->toString());
            break;
        case 2:
            this->setupTableData(var);
            break;
        case 3:
            this->setupPiecewiseData(var);
            break;
        default:
            qWarning("%s", QString(tr("Dialog: Cannot display data of an unknown type!")).toStdString().c_str());
    }
}

void VariableDialog::setupTableData(const Variable& var)
{
    TableType* table = (TableType*)var.data();
    for (auto it = table->data().cbegin(); it != table->data().cend(); ++it) {
        QList<QStandardItem*> list;
        QStandardItem* cell1 = new QStandardItem(it->first);
        QStandardItem* cell2 = new QStandardItem(it->second);
        list << cell1 << cell2;
        this->tableModel->appendRow(list);
    }
}

DataType* VariableDialog::collectTableData() const
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

void VariableDialog::setupPiecewiseData(const Variable& var)
{
    PiecewiseFunctionType* table = (PiecewiseFunctionType*)var.data();
    for (auto it = table->data().cbegin(); it != table->data().cend(); ++it) {
        QList<QStandardItem*> list;
        QStandardItem* cell1 = new QStandardItem(it->at(0));
        QStandardItem* cell2 = new QStandardItem(it->at(1));
        QStandardItem* cell3 = new QStandardItem(it->at(2));
        list << cell1 << cell2 << cell3;
        this->piecewiseModel->appendRow(list);
    }
}

DataType* VariableDialog::collectPiecewiseData() const
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

void VariableDialog::on_btnTableAdd_pressed()
{
    QList<QStandardItem*> list;
    QStandardItem* cell1 = new QStandardItem("0");
    QStandardItem* cell2 = new QStandardItem("0");
    list << cell1 << cell2;
    this->tableModel->appendRow(list);
}

void VariableDialog::on_btnTableDel_pressed()
{
    QModelIndexList indices = ui->tblTable->selectionModel()->selectedIndexes();

    foreach (QModelIndex index, indices) {
        this->tableModel->removeRow(index.row());
    }
}

void VariableDialog::on_btnPiecewiseAdd_pressed()
{
    QList<QStandardItem*> list;
    QStandardItem* cell1 = new QStandardItem("0");
    QStandardItem* cell2 = new QStandardItem("1");
    QStandardItem* cell3 = new QStandardItem("x");
    list << cell1 << cell2 << cell3;
    this->piecewiseModel->appendRow(list);
}

void VariableDialog::on_btnPiecewiseDel_pressed()
{
    QModelIndexList indices = ui->tblPiecewise->selectionModel()->selectedIndexes();

    foreach (QModelIndex index, indices) {
        this->piecewiseModel->removeRow(index.row());
    }
}

void VariableDialog::visit(ConstantType& type)
{
    this->ui->cmbType->setCurrentIndex(0);
}

void VariableDialog::visit(FunctionType& type)
{
    this->ui->cmbType->setCurrentIndex(1);
}

void VariableDialog::visit(TableType& type)
{
    this->ui->cmbType->setCurrentIndex(2);
}

void VariableDialog::visit(PiecewiseFunctionType& type)
{
    this->ui->cmbType->setCurrentIndex(3);
}

void VariableDialog::visit(VariableLinkType& type)
{
}
