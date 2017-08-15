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

    QRegExpValidator* varNameValidator = new QRegExpValidator(QRegExp(REGEXPR_VARNAME), this);
    ui->editName->setValidator(varNameValidator);

    this->varDict = varDict;

    this->dataWidget = new DataTypeWidget(varDict, ui->frData);
    ui->frData->layout()->addWidget(this->dataWidget);
}

VariableDialog::VariableDialog(const Variable& var, const QHash<QString, Variable>& varDict, QWidget *parent) :
    VariableDialog(varDict, parent)
{
    this->varDict.remove(var.name());
    ui->editName->setText(var.name());
    this->dataWidget->setDataType(var.data());
}

VariableDialog::~VariableDialog()
{
    delete ui;
}

Variable VariableDialog::data()
{
    QString name = ui->editName->text();
    DataType* data = this->dataWidget->data();

    return Variable(name, data);
}

void VariableDialog::accept()
{
    if (this->dataWidget->check() == -1)
        return;
    QDialog::accept();
}
