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
        QAbstractItemModel* model = ui->listVariables->model();
        int row = model->rowCount();
        model->insertRows(row, 1);
        QModelIndex index = model->index(row, 0);
        QVariant data = dialog->getData();
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
        QVariant result = dialog->getData();
        ui->listVariables->model()->setData(index, result);
    }
}
