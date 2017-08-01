#ifndef VARIABLEDIALOG_H
#define VARIABLEDIALOG_H

#include <QDialog>
#include <QStandardItemModel>

#include "../data/variable.h"

namespace Ui {
class VariableDialog;
}

class VariableDialog : public QDialog
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

private:
    Ui::VariableDialog *ui;
    QHash<QString, Variable> varDict;
    QStandardItemModel* tableModel;

    void setData(const Variable& var);
    DataType* collectTableData() const;
    void setupTableData(const Variable& var);
};

#endif // VARIABLEDIALOG_H
