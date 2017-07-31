#include "variabledialog.h"
#include "ui_variabledialog.h"

#include "../../data/common.h"
#include <QMessageBox>

VariableDialog::VariableDialog(const QHash<QString, Variable>& varDict, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VariableDialog)
{
    ui->setupUi(this);

    ui->cmbType->addItem(tr("Constant"));
    ui->cmbType->addItem(tr("Function"));
    ui->cmbType->addItem(tr("Table"));
    ui->cmbType->addItem(tr("Piecewise function"));
    ui->cmbType->setCurrentIndex(0);

    QRegExpValidator* varNameValidator = new QRegExpValidator(QRegExp(REGEXPR_VARNAME), this);
    ui->editName->setValidator(varNameValidator);

    QRegExpValidator* constantValidator = new QRegExpValidator(QRegExp(REGEXPR_DOUBLE), this);
    ui->editConstant->setValidator(constantValidator);

    this->varDict = varDict;
}

VariableDialog::VariableDialog(const Variable& var, const QHash<QString, Variable>& varDict, QWidget *parent) :
    VariableDialog(varDict, parent)
{
    ui->editName->setText(var.name());
    ui->editConstant->setText(var.data()->toString());
    ui->cmbType->setCurrentIndex(var.data()->type());
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
        case DTLib::CONSTANT:
            data = new ConstantType(ui->editConstant->text());
            break;
        default:
            qWarning("%s", QString(tr("VariableDialog: Unknown type!")).toStdString().c_str());
    }

    return Variable(name, data);
}

void VariableDialog::on_cmbType_currentIndexChanged(int index)
{
    switch(index)
    {
        case DTLib::CONSTANT:
            ui->editConstant->clear();
            ui->editConstant->show();
            ui->editFunction->hide();
            break;
        case DTLib::FUNCTION:
            ui->editConstant->hide();
            ui->editFunction->clear();
            ui->editFunction->show();
            break;
        default:
            qWarning("%s", QString(tr("VariableDialog: Unknown type!")).toStdString().c_str());
    }
}

void VariableDialog::accept()
{
    if (this->varDict.contains(ui->editName->text()))
    {
        QMessageBox::warning(this, tr("Error"), tr("Chosen variable name is already in use!"));
    }
    else
    {
        QDialog::accept();
    }
}
