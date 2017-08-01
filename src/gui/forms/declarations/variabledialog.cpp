#include "variabledialog.h"
#include "ui_variabledialog.h"

#include "../../data/common.h"
#include <QMessageBox>
#include "../../expression.h"

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
    this->varDict.remove(var.name());
    ui->editName->setText(var.name());
    ui->cmbType->setCurrentIndex(var.data()->type());
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
        case DTLib::CONSTANT:
            data = new ConstantType(ui->editConstant->text());
            break;
        case DTLib::FUNCTION:
            data = new FunctionType(ui->editFunction->text());
            break;
        default:
            qWarning("%s", QString(tr("VariableDialog: Cannot retrieve data from an unknown type!")).toStdString().c_str());
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
            qWarning("%s (%d)", QString(tr("VariableDialog: Unknown type selected!")).toStdString().c_str(), index);
    }
}

void VariableDialog::accept()
{
    if (this->varDict.contains(ui->editName->text()))
    {
        QMessageBox::warning(this, tr("Error"), tr("Chosen variable name is already in use!"));
        return;
    }
    else if (ui->cmbType->currentIndex() == DTLib::FUNCTION)
    {
        if (ui->editFunction->text().isEmpty())
        {
            QMessageBox::warning(this, tr("Error"), tr("Function formula not entered!"));
            return;
        }
        std::vector<std::string> vars;
        vars.push_back("x");
        vars.push_back("y");
        vars.push_back("z");
        vars.push_back("time");
        vars.push_back("temperature");
        if (!Expression::isValid(ui->editFunction->text().toStdString(), vars))
        {
            QMessageBox::warning(this, tr("Error"), tr("Incorrect format of function formula!"));
            return;
        }
    }

    QDialog::accept();
}

void VariableDialog::setData(const Variable &var)
{
    switch(ui->cmbType->currentIndex())
    {
        case DTLib::CONSTANT:
            ui->editConstant->setText(var.data()->toString());
            break;
        case DTLib::FUNCTION:
            ui->editFunction->setText(var.data()->toString());
            break;
        default:
            qWarning("%s", QString(tr("VariableDialog: Cannot display data of an unknown type!")).toStdString().c_str());
    }
}
