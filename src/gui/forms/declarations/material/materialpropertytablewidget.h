#ifndef MATERIALPROPERTYTABLEWIDGET_H
#define MATERIALPROPERTYTABLEWIDGET_H

#include <QWidget>
#include <QVector>
#include "../../../data/datatype.h"
#include "../../../data/material/materialproperty.h"

namespace Ui {
class MaterialPropertyTableWidget;
}

class MaterialPropertyTableWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MaterialPropertyTableWidget(QWidget *parent = 0, bool withHeader = true);
    ~MaterialPropertyTableWidget();

    void addProperty(const MaterialProperty* property);
    void addRow(const QString& name, DataType* data, const QString& unit, const QString& symbol);

private:
    Ui::MaterialPropertyTableWidget *ui;

    void createHeader();
};

#endif // MATERIALPROPERTYTABLEWIDGET_H
