#ifndef DECLARATIONSWIDGET_H
#define DECLARATIONSWIDGET_H

#include <QWidget>
#include <QMenu>
#include <QStandardItemModel>
#include <QString>
#include "../data/variable.h"

namespace espreso
{
    namespace Ui {
    class DeclarationsWidget;
    }

    class DeclarationsWidget : public QWidget
    {
        Q_OBJECT

    public:
        explicit DeclarationsWidget(QWidget *parent = 0);
        ~DeclarationsWidget();

    private slots:
        void on_DeclarationTree_customContextMenuRequested(const QPoint &pos);
        void treeNewItem();
        void treeEditItem();
        void treeDelItem();

        void on_DeclarationTree_doubleClicked(const QModelIndex &index);

    private:
        Ui::DeclarationsWidget *ui;

        QAction* m_newItem;
        QAction* m_editItem;
        QAction* m_delItem;

        int m_varRow = -1;
        int m_csRow = -1;
        int m_matRow = -1;

        QStandardItemModel* m_treeModel;
        QStandardItem* m_treeNodeVars;
        QStandardItem* m_treeNodeCS;
        QStandardItem* m_treeNodeMats;

        QVector<Variable> m_variables;
        QHash<QString, Variable> m_varDict;

        void setupTree();
        void createActions();
        void createEditDialog(const QModelIndex& item);
    };
}
#endif // DECLARATIONSWIDGET_H
