#include "variabledialog.h"
#include "ui_variabledialog.h"

VariableDialog::VariableDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VariableDialog)
{
    ui->setupUi(this);
}

VariableDialog::~VariableDialog()
{
    delete ui;
}

Variable VariableDialog::getVariable()
{
    QString name = ui->editName->text();
    StringType* data = new StringType(ui->editData->text());

    return Variable(name, data);
}
