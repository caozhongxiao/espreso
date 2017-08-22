#ifndef MATERIALPROPERTYTABLEWIDGET_H
#define MATERIALPROPERTYTABLEWIDGET_H

#include <QWidget>
#include "../../../data/datatype.h"

namespace Ui {
class MaterialPropertyTableWidget;
}

class MaterialPropertyTableWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MaterialPropertyTableWidget(QWidget *parent = 0);
    ~MaterialPropertyTableWidget();

    void addRow(const QString& name, DataType* data, const QString& unit, const QString& abbrev);

private:
    Ui::MaterialPropertyTableWidget *ui;

};

#endif // MATERIALPROPERTYTABLEWIDGET_H
