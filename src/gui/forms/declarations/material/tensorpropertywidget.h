#ifndef TENSORPROPERTYWIDGET_H
#define TENSORPROPERTYWIDGET_H

#include <QWidget>

namespace Ui {
class TensorPropertyWidget;
}

class TensorPropertyWidget : public QWidget
{
    Q_OBJECT

public:
    explicit TensorPropertyWidget(QWidget *parent = 0);
    ~TensorPropertyWidget();

private:
    Ui::TensorPropertyWidget *ui;
};

#endif // TENSORPROPERTYWIDGET_H
