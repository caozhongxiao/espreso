#ifndef VARIABLEDIALOG_H
#define VARIABLEDIALOG_H

#include <QDialog>
#include <QStandardItemModel>

#include "../data/datatype.h"
#include "../data/datatype.h"
#include "datatypewidget.h"

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

private:
    Ui::VariableDialog *ui;
    QHash<QString, Variable> varDict;

    DataTypeWidget* dataWidget;

    void setData(const Variable& var);
};

#endif // VARIABLEDIALOG_H
