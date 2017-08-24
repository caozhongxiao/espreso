#include "tablewidget.h"
#include "ui_tablewidget.h"

TableWidget::TableWidget(int columns, const QStringList& headlines, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TableWidget)
{
    ui->setupUi(this);

    this->mTable = this->ui->table;

    this->mModel = new QStandardItemModel(this);
    ui->table->setModel(mModel);

    this->mModel->setColumnCount(columns);

    this->mModel->setHorizontalHeaderLabels(headlines);

    this->addCleanRow();
}

TableWidget::~TableWidget()
{
    delete ui;
}

void TableWidget::addCleanRow()
{
    QList<QStandardItem*> row;
    for (int i = 0; i < mModel->rowCount(); ++i)
    {
        row << new QStandardItem();
    }
    this->mModel->appendRow(row);
}

void TableWidget::addRow(const QList<QString>& rowData)
{
    QList<QStandardItem*> row;
    foreach (QString data, rowData) {
        row.append(new QStandardItem(data));
    }
    this->mModel->insertRow(mModel->rowCount() - 1, row);
}

void TableWidget::addData(const QList<QList<QString> >& data)
{
    foreach (QList<QString> row, data) {
        this->addRow(row);
    }
}
