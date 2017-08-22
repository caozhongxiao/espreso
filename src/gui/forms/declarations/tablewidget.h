#ifndef TABLEWIDGET_H
#define TABLEWIDGET_H

#include <QWidget>
#include <QStandardItemModel>
#include <QStringList>
#include "../validators/validatorfactory.h"

namespace Ui {
class TableWidget;
}

class TableWidget : public QWidget
{
    Q_OBJECT

public:
    explicit TableWidget(int columns, const QStringList& headlines,
                         const QVector<ValidatorFactory*>& validators, QWidget *parent = 0);
    ~TableWidget();

    void addRow(const QList<QString>& rowData);
    void addData(const QList<QList<QString> >& data);

private:
    Ui::TableWidget *ui;

    QStandardItemModel* mModel;
};

#endif // TABLEWIDGET_H
