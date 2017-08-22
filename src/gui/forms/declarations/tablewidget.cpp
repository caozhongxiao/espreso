#include "tablewidget.h"
#include "ui_tablewidget.h"

#include "../validators/validatordelegate.h"

TableWidget::TableWidget(int columns, const QStringList& headlines,
                         const QVector<ValidatorFactory*>& validators, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TableWidget)
{
    ui->setupUi(this);

    this->mModel = new QStandardItemModel(this);
    ui->table->setModel(mModel);

    this->mModel->setColumnCount(columns);

    int col = 0;
    foreach (ValidatorFactory* vf, validators) {
        ui->table->setItemDelegateForColumn(col, new ValidatorDelegate(vf, ui->table));
        col++;
    }

    this->mModel->setHorizontalHeaderLabels(headlines);
}

TableWidget::~TableWidget()
{
    delete ui;
}

void TableWidget::addRow(const QList<QString>& rowData)
{
    QList<QStandardItem*> row;
    foreach (QString data, rowData) {
        row.append(new QStandardItem(data));
    }
    this->mModel->appendRow(row);
}

void TableWidget::addData(const QList<QList<QString> >& data)
{
    foreach (QList<QString> row, data) {
        this->addRow(row);
    }
}
