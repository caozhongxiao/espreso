#ifndef VARIABLEDIALOG_H
#define VARIABLEDIALOG_H

#include <QDialog>
#include <QStandardItemModel>

#include "../data/variable.h"
#include "../data/datatype.h"

namespace Ui {
class VariableDialog;
}

class VariableDialog : public QDialog, public DataTypeVisitor
{
    Q_OBJECT

public:
    explicit VariableDialog(const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    VariableDialog(const Variable& var, const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    //VariableDialog(QVariant data, QWidget *parent = 0);
    ~VariableDialog();
    Variable data();
    void accept() override;

private slots:
    void on_cmbType_currentIndexChanged(int index);

    void on_btnTableAdd_pressed();

    void on_btnTableDel_pressed();

    void on_btnPiecewiseAdd_pressed();

    void on_btnPiecewiseDel_pressed();

private:
    Ui::VariableDialog *ui;
    QHash<QString, Variable> varDict;
    QStandardItemModel* tableModel;
    QStandardItemModel* piecewiseModel;
    std::vector<std::string> fnVars;

    void setData(const Variable& var);
    DataType* collectTableData() const;
    void setupTableData(const Variable& var);
    DataType* collectPiecewiseData() const;
    void setupPiecewiseData(const Variable& var);

    void visit(const ConstantType& type) override;
    void visit(const FunctionType& type) override;
    void visit(const TableType& type) override;
    void visit(const PiecewiseFunctionType& type) override;
    void visit(const VariableLinkType& type) override;
};

#endif // VARIABLEDIALOG_H
