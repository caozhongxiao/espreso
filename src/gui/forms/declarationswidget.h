#ifndef DECLARATIONSWIDGET_H
#define DECLARATIONSWIDGET_H

#include <QWidget>
#include <QMenu>
#include <QStandardItemModel>

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

private:
    Ui::DeclarationsWidget *ui;

    QAction* newItem;
    QAction* editItem;
    QAction* delItem;

    qint32 varCount = 0;
    qint32 csCount = 0;
    qint32 matCount = 0;

    QStandardItemModel* treeModel;
    QStandardItem* treeVars;

    void setupTree();
    void createActions();
};

#endif // DECLARATIONSWIDGET_H
