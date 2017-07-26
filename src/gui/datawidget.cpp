#include "datawidget.h"
#include "ui_datawidget.h"

#include "variablemodel.h"
#include "variabledialog.h"

DataWidget::DataWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DataWidget)
{
    ui->setupUi(this);

    // VARIABLES
    ui->listVariables->setModel(new VariableListModel(this));
}

DataWidget::~DataWidget()
{
    delete ui;
}

void DataWidget::on_btnVarAdd_pressed()
{
    VariableDialog* dialog = new VariableDialog(this);
    if (dialog->exec() == QDialog::Accepted)
    {
        Variable var = dialog->getVariable();
        QAbstractItemModel* model = ui->listVariables->model();
        int row = model->rowCount();
        model->insertRows(row, 1);
        QModelIndex index = model->index(row, 0);
        QVariant data;
        data.setValue(var);
        model->setData(index, data);
    }
}
