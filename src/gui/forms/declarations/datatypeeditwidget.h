#ifndef DATATYPEEDITWIDGET_H
#define DATATYPEEDITWIDGET_H

#include <QWidget>

#include <QLineEdit>
#include <QComboBox>
#include "tabletypewidget.h"
#include "piecewisetypewidget.h"
#include "../../data/datatype.h"
#include "../elements/expressionedit.h"

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
    bool isValid();

private slots:
    void changeType(int index);

private:
    Ui::DataTypeEditWidget *ui;

    ExpressionEdit* uiExpression;
    TableTypeWidget* uiTable;
    PiecewiseTypeWidget* uiPiecewise;

    int activeType;

    void createUi();
    void initExpression(const ExpressionType*);
    void initTable(const TableType*);
    void initPiecewise(const PiecewiseFunctionType*);
};

#endif // DATATYPEEDITWIDGET_H
