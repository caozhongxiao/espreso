#ifndef VARIABLEDIALOG_H
#define VARIABLEDIALOG_H

#include <QDialog>

#include "variablemodel.h"

namespace Ui {
class VariableDialog;
}

class VariableDialog : public QDialog
{
    Q_OBJECT

public:
    explicit VariableDialog(QWidget *parent = 0);
    ~VariableDialog();
    Variable getVariable();

private:
    Ui::VariableDialog *ui;
};

#endif // VARIABLEDIALOG_H
