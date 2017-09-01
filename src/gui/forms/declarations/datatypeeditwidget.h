#ifndef DATATYPEEDITWIDGET_H
#define DATATYPEEDITWIDGET_H

#include <QWidget>

#include <QLineEdit>
#include <QComboBox>
#include "tabletypewidget.h"
#include "piecewisetypewidget.h"
#include "../../../config/configuration.h"
#include "../elements/expressionedit.h"

using namespace espreso;

namespace Ui {
class DataTypeEditWidget;
}

class DataTypeEditWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DataTypeEditWidget(QWidget *parent = 0);
    DataTypeEditWidget(const ECFValue& data, QWidget* parent = 0);
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
    void initExpression(const ECFValue&);
    void initTable(const ECFValue&);
    void initPiecewise(const ECFValue&);
};

#endif // DATATYPEEDITWIDGET_H
