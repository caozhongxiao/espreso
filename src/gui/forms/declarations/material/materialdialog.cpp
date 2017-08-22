#include "materialdialog.h"
#include "ui_materialdialog.h"

MaterialDialog::MaterialDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MaterialDialog)
{
    ui->setupUi(this);
}

MaterialDialog::~MaterialDialog()
{
    delete ui;
}
