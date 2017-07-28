#ifndef DATAWIDGET_H
#define DATAWIDGET_H

#include <QWidget>

namespace Ui {
class DataWidget;
}

class DataWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DataWidget(QWidget *parent = 0);
    ~DataWidget();

private slots:
    void on_btnVarAdd_pressed();

    void on_btnVarDel_pressed();

    void on_listVariables_doubleClicked(const QModelIndex &index);

    void on_chbVariables_stateChanged(int arg1);

private:
    Ui::DataWidget *ui;
};

#endif // DATAWIDGET_H
