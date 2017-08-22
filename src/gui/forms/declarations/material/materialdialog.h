#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include <QDialog>

namespace Ui {
class MaterialDialog;
}

class MaterialDialog : public QDialog
{
    Q_OBJECT

public:
    explicit MaterialDialog(QWidget *parent = 0);
    ~MaterialDialog();

private:
    Ui::MaterialDialog *ui;
};

#endif // MATERIALDIALOG_H
