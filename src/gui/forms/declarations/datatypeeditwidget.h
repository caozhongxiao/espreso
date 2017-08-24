#ifndef DATATYPEEDITWIDGET_H
#define DATATYPEEDITWIDGET_H

#include <QWidget>

#include <QLineEdit>
#include <QComboBox>
#include "tabletypewidget.h"
#include "piecewisetypewidget.h"
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

signals:
    void validStateChanged(bool valid);

private slots:
    void changeType(int index);
    void changeValidState(bool valid);

private:
    Ui::DataTypeEditWidget *ui;

    QLineEdit* uiExpression;
    TableTypeWidget* uiTable;
    PiecewiseTypeWidget* uiPiecewise;

    int activeType;

    void createUi();
    void initExpression(const ExpressionType*);
    void initTable(const TableType*);
    void initPiecewise(const PiecewiseFunctionType*);
};

#endif // DATATYPEEDITWIDGET_H
