#ifndef MATRIXPROPERTYWIDGET_H
#define MATRIXPROPERTYWIDGET_H

#include <QWidget>
#include <QHash>
#include <QVector>
#include <QStandardItemModel>

#include "../data/variable.h"
#include "../data/datatype.h"
#include "../data/material.h"
#include "datatypewidget.h"

namespace Ui {
class MatrixPropertyWidget;
}

class MatrixPropertyWidget : public QWidget, public MaterialPropertyVisitor
{
    Q_OBJECT

public:
    explicit MatrixPropertyWidget(const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    MatrixPropertyWidget(MaterialProperty* property, const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    ~MatrixPropertyWidget();

    void setData(MaterialProperty* property);
    MaterialProperty* data();

    void visit(BasicProperty& p) override;
    void visit(IsotropicProperty& p) override;
    void visit(DiagonalProperty& p) override;
    void visit(SymmetricProperty& p) override;
    void visit(AnisotropicProperty& p) override;

private slots:
    void on_btnSave_pressed();

private:
    Ui::MatrixPropertyWidget *ui;

    QHash<QString, Variable> varDict;
    QStandardItemModel* matrixModel;
    int propertyType;
    MaterialProperty* mProperty;
    DataTypeWidget* matrixItemWidget;
    QVector<QVector<DataType*> > matrix;

    int pRow;
    int pColumn;
    void matrixBtnPressed(int row, int column);
};

#endif // MATRIXPROPERTYWIDGET_H
