#ifndef DECLARATIONSWIDGET_H
#define DECLARATIONSWIDGET_H

#include <QWidget>

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


    void setupTree();
    void createActions();
};

#endif // DECLARATIONSWIDGET_H
