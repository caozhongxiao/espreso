#include "variabledialog.h"
#include "ui_variabledialog.h"

VariableDialog::VariableDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VariableDialog)
{
    ui->setupUi(this);
}

VariableDialog::VariableDialog(QVariant data, QWidget *parent) :
    VariableDialog(parent)
{
    if (!data.canConvert<Variable>()) return;
    Variable var = data.value<Variable>();
    ui->editName->setText(var.getName());
    ui->editData->setText(var.getData()->toString());
}

VariableDialog::~VariableDialog()
{
    delete ui;
}

QVariant VariableDialog::getData()
{
    QString name = ui->editName->text();
    StringType* data = new StringType(ui->editData->text());
    Variable v(name, data);
    QVariant output;
    output.setValue(v);

    return output;
}
