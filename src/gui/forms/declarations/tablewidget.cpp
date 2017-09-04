#include "tablewidget.h"
#include "ui_tablewidget.h"

#include <QtDebug>

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

    connect(mTable, &QTableView::doubleClicked,
            this, &TableWidget::onItemDoubleClick);

    this->mTable->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(mTable, &QTableView::customContextMenuRequested,
            this, &TableWidget::onContextMenu);

    this->mActionDelete = new QAction(tr("&Delete"), this);
    connect(mActionDelete, &QAction::triggered, this, &TableWidget::deleteItem);
}

TableWidget::~TableWidget()
{
    delete ui;
}

void TableWidget::addCleanRow()
{
    QList<QStandardItem*> row;
    for (int i = 0; i < mModel->columnCount(); ++i)
    {
        //row << new QStandardItem(columnDefaultValue(i));
        row << new QStandardItem();
    }
    this->mModel->appendRow(row);
}

void TableWidget::addRow(const QVector<QString>& rowData)
{
    QList<QStandardItem*> row;
    foreach (QString data, rowData) {
        row.append(new QStandardItem(data));
    }
    this->mModel->insertRow(mModel->rowCount() - 1, row);
}

void TableWidget::addData(const QVector<QVector<QString> >& data)
{
    foreach (QVector<QString> row, data) {
        this->addRow(row);
    }
}

void TableWidget::onItemDoubleClick(const QModelIndex &index)
{
    if (index.row() != this->mModel->rowCount() - 1)
        return;

    this->addCleanRow();
}

void TableWidget::onContextMenu(const QPoint& point)
{
    this->toDelete = this->mTable->indexAt(point);
    QMenu menu(this);
    menu.addAction(this->mActionDelete);
    menu.exec(ui->table->mapToGlobal(point));
}

void TableWidget::deleteItem()
{
    if (toDelete.row() == -1 || toDelete.column() == -1)
        return;
    if (toDelete.row() == mModel->rowCount() - 1)
        return;

    this->mModel->removeRow(toDelete.row());
}

bool TableWidget::isValid()
{
    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        for (int col = 0; col < mModel->columnCount(); ++col)
        {
            QModelIndex index = this->mModel->index(row, col);
            QString val = index.data().toString();
            if (val.isEmpty())
                return false;
        }
    }

    return true;
}
