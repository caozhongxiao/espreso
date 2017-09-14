#ifndef DATATYPEWIDGET_H
#define DATATYPEWIDGET_H

#include <QWidget>
#include <QStandardItemModel>

#include "../../data/datatype.h"
#include "../../data/variable.h"
#include "../../data/common.h"

namespace Ui {
class DataTypeWidget;
}

class DataTypeWidget : public QWidget, public DataTypeVisitor
{
    Q_OBJECT

public:
    explicit DataTypeWidget(const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    DataTypeWidget(DataType* data, const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    ~DataTypeWidget();
    DataType* data() const;
    int check();
    DataType* checkedData();
    void setDataType(DataType*);

    virtual void visit(ConstantType& type) override;
    virtual void visit(FunctionType& type) override;
    virtual void visit(ExpressionType& type) override {}
    virtual void visit(TableType& type) override;
    virtual void visit(PiecewiseFunctionType& type) override;
    virtual void visit(VariableLinkType& type) override;
    virtual void visit(DummyType& type) override {}

private slots:
    void on_cmbType_currentIndexChanged(int index);

    void on_btnTableAdd_pressed();

    void on_btnTableDel_pressed();

    void on_btnPiecewiseAdd_pressed();

    void on_btnPiecewiseDel_pressed();

private:
    Ui::DataTypeWidget *ui;

    QHash<QString, Variable> varDict;
    QStandardItemModel* tableModel;
    QStandardItemModel* piecewiseModel;

    void setData(DataType* data);
    DataType* collectTableData() const;
    void setupTableData(DataType* data);
    DataType* collectPiecewiseData() const;
    void setupPiecewiseData(DataType* data);
    void setupVariableLinkData(DataType* data);

protected:
    virtual void setupCheckbox();
};

#endif // DATATYPEWIDGET_H
