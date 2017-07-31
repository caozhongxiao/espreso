#include "variabledialog.h"
#include "ui_variabledialog.h"

VariableDialog::VariableDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VariableDialog)
{
    ui->setupUi(this);
}

VariableDialog::VariableDialog(Variable var, QWidget *parent) :
    VariableDialog(parent)
{
    ui->editName->setText(var.name());
    ui->editData->setText(var.data()->toString());
}

VariableDialog::VariableDialog(QVariant data, QWidget *parent) :
    VariableDialog(parent)
{
    if (!data.canConvert<Variable>())
    {
        qWarning("%s", tr("VariableDialog: No data provided for the dialog!").toStdString().c_str());
        return;
    }
    Variable var = data.value<Variable>();
    ui->editName->setText(var.name());
    ui->editData->setText(var.data()->toString());
}

VariableDialog::~VariableDialog()
{
    delete ui;
}

Variable VariableDialog::data()
{
    QString name = ui->editName->text();
    StringType* data = new StringType(ui->editData->text());

    return Variable(name, data);
}
