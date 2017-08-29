#ifndef MATERIALPROPERTIESWIDGET_H
#define MATERIALPROPERTIESWIDGET_H

#include <QWidget>
#include <QVector>

#include "tensorpropertywidget.h"
#include "materialpropertytablewidget.h"
#include "../../../data/material/scalarproperty.h"

namespace Ui {
class MaterialPropertiesWidget;
}

class MaterialPropertiesWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MaterialPropertiesWidget(const QVector<TensorProperty>& tensors,
                                      const QVector<ScalarProperty>& scalars,
                                      QWidget *parent = 0);
    ~MaterialPropertiesWidget();

private:
    Ui::MaterialPropertiesWidget *ui;

    QVector<TensorPropertyWidget*> tensorWidgets;
    MaterialPropertyTableWidget* scalarWidget;

};

#endif // MATERIALPROPERTIESWIDGET_H
