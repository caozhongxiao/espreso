#ifndef DATATYPEEDITWIDGET_H
#define DATATYPEEDITWIDGET_H

#include <QWidget>

#include <QLineEdit>
#include <QComboBox>
#include "tablewidget.h"
#include "../../data/datatype.h"

namespace Ui {
class DataTypeEditWidget;
}

class DataTypeEditWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DataTypeEditWidget(QWidget *parent = 0);
    DataTypeEditWidget(const DataType* data, QWidget* parent = 0);
    ~DataTypeEditWidget();

    QComboBox* createComboBox(QWidget* parent = nullptr);

private slots:
    void typeChanged(int index);

private:
    Ui::DataTypeEditWidget *ui;

    QLineEdit* uiExpression;
    TableWidget* uiTable;
    TableWidget* uiPiecewise;

    int activeType;

    void createUi();
    void initExpression(const ExpressionType*);
    void initTable(const TableType*);
    void initPiecewise(const PiecewiseFunctionType*);
};

#endif // DATATYPEEDITWIDGET_H
