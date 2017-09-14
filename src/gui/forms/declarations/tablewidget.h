#ifndef TABLEWIDGET_H
#define TABLEWIDGET_H

#include <QWidget>
#include <QStandardItemModel>
#include <QStringList>
#include <QTableView>
#include "../validators/validatorfactory.h"
#include <QAction>
#include <QMenu>

namespace espreso
{

    namespace Ui {
    class TableWidget;
    }

    class TableWidget : public QWidget
    {
        Q_OBJECT

    public:
        virtual ~TableWidget();

        virtual void addRow(const QVector<QString>& rowData);
        virtual void addData(const QVector<QVector<QString> >& data);
        virtual void addData(const QString& data) = 0;
        virtual bool isValid();

    protected:
        QTableView* mTable;
        QStandardItemModel* mModel;

        explicit TableWidget(int columns, const QStringList& headlines, QWidget *parent = 0);

        virtual QString columnDefaultValue(int) const { return ""; }

    private slots:
        void onItemDoubleClick(const QModelIndex& index);
        void onContextMenu(const QPoint& point);
        void deleteItem();

    private:
        Ui::TableWidget *ui;

        QModelIndex toDelete;
        QAction* mActionDelete;

        void addCleanRow();
    };
}
#endif // TABLEWIDGET_H
