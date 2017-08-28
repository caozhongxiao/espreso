#ifndef TENSORPROPERTYWIDGET_H
#define TENSORPROPERTYWIDGET_H

#include <QWidget>
#include "../../../data/material/tensorproperty.h"
#include "../../../data/material/tensorpropertymodel.h"
#include "materialpropertytablewidget.h"

namespace Ui {
class TensorPropertyWidget;
}

class TensorPropertyWidget : public QWidget
{
    Q_OBJECT

public:
    explicit TensorPropertyWidget(const TensorProperty& property, QWidget *parent = 0);
    ~TensorPropertyWidget();

private slots:
    void onIndexChanged(int index);

private:
    Ui::TensorPropertyWidget *ui;

    TensorProperty mProperty;
    QVector<MaterialPropertyTableWidget*> mWidgets;
};

#endif // TENSORPROPERTYWIDGET_H
