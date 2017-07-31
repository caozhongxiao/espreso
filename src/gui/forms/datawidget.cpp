#include "datawidget.h"
#include "ui_datawidget.h"

#include "declarations/variabledialog.h"
#include "../models/variablemodel.h"

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
        QAbstractItemModel* model = ui->listVariables->model();
        int row = model->rowCount();
        model->insertRows(row, 1);
        QModelIndex index = model->index(row, 0);
        Variable var = dialog->data();
        QVariant data;
        data.setValue(var);
        model->setData(index, data);
    }
}

void DataWidget::on_btnVarDel_pressed()
{
    QModelIndexList selected = ui->listVariables->selectionModel()->selectedIndexes();
    QAbstractItemModel* model = ui->listVariables->model();
    for (int i = selected.count() - 1; i > -1; --i)
    {
        model->removeRow(selected.at(i).row());
    }
}

void DataWidget::on_listVariables_doubleClicked(const QModelIndex &index)
{
    QVariant data = index.data(Qt::EditRole);
    VariableDialog* dialog = new VariableDialog(data, this);
    if (dialog->exec() == QDialog::Accepted)
    {
        Variable var = dialog->data();
        QVariant result;
        result.setValue(var);
        ui->listVariables->model()->setData(index, result);
    }
}

void DataWidget::on_chbVariables_stateChanged(int arg1)
{
    if (ui->chbVariables->isChecked())
    {
        ui->listVariables->hide();
        ui->btnVarAdd->hide();
        ui->btnVarDel->hide();
    }
    else
    {
        ui->listVariables->show();
        ui->btnVarAdd->show();
        ui->btnVarDel->show();
    }
}
