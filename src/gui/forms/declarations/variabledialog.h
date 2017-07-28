#ifndef VARIABLEDIALOG_H
#define VARIABLEDIALOG_H

#include <QDialog>

#include "../data/variable.h"

namespace Ui {
class VariableDialog;
}

class VariableDialog : public QDialog
{
    Q_OBJECT

public:
    explicit VariableDialog(QWidget *parent = 0);
    VariableDialog(QVariant data, QWidget *parent = 0);
    ~VariableDialog();
    QVariant data();

private:
    Ui::VariableDialog *ui;
};

#endif // VARIABLEDIALOG_H
