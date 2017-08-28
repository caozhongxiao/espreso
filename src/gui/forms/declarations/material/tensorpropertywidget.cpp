#include "tensorpropertywidget.h"
#include "ui_tensorpropertywidget.h"

TensorPropertyWidget::TensorPropertyWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TensorPropertyWidget)
{
    ui->setupUi(this);
}

TensorPropertyWidget::~TensorPropertyWidget()
{
    delete ui;
}
