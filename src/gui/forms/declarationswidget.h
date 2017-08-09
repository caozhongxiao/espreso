#ifndef DECLARATIONSWIDGET_H
#define DECLARATIONSWIDGET_H

#include <QWidget>
#include <QMenu>
#include <QStandardItemModel>
#include <QString>
#include "../data/variable.h"
#include "../data/material.h"

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

    QAction* newItem;
    QAction* editItem;
    QAction* delItem;

    QStandardItemModel* treeModel;
    QStandardItem* treeNodeVars;
    QStandardItem* treeNodeCS;
    QStandardItem* treeNodeMats;

    QVector<Variable> variables;
    QHash<QString, Variable> varDict;

    QVector<Material> materials;
    QHash<QString, Material> matDict;

    void setupTree();
    void createActions();
    void createEditDialog(const QModelIndex& item);
};

#endif // DECLARATIONSWIDGET_H
