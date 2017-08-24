#ifndef TABLEWIDGET_H
#define TABLEWIDGET_H

#include <QWidget>
#include <QStandardItemModel>
#include <QStringList>
#include <QTableView>
#include "../validators/validatorfactory.h"

namespace Ui {
class TableWidget;
}

class TableWidget : public QWidget
{
    Q_OBJECT

public:
    explicit TableWidget(int columns, const QStringList& headlines, QWidget *parent = 0);
    virtual ~TableWidget();

    virtual void addRow(const QList<QString>& rowData);
    virtual void addData(const QList<QList<QString> >& data);

protected:
    QTableView* mTable;
    QStandardItemModel* mModel;

private slots:
    void onItemDoubleClick(const QModelIndex& index);

private:
    Ui::TableWidget *ui;

    void addCleanRow();
};

#endif // TABLEWIDGET_H
