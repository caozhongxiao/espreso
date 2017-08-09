#ifndef DATATYPEWIDGET_H
#define DATATYPEWIDGET_H

#include <QWidget>
#include <QStandardItemModel>

#include "../data/datatype.h"
#include "../data/variable.h"
#include "../data/common.h"

namespace Ui {
class DataTypeWidget;
}

class DataTypeWidget : public QWidget, public DataTypeVisitor
{
    Q_OBJECT

public:
    explicit DataTypeWidget(const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    ~DataTypeWidget();
    DataType* data() const;
    int check();

    virtual void visit(const ConstantType& type) override;
    virtual void visit(const FunctionType& type) override;
    virtual void visit(const TableType& type) override;
    virtual void visit(const PiecewiseFunctionType& type) override;
    virtual void visit(const VariableLinkType& type) override;

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

    void setData(const DataType* data);
    DataType* collectTableData() const;
    void setupTableData(const DataType* data);
    DataType* collectPiecewiseData() const;
    void setupPiecewiseData(const DataType* data);
    void setupVariableLinkData(const DataType* data);

protected:
    virtual void setupCheckbox();
};

#endif // DATATYPEWIDGET_H
