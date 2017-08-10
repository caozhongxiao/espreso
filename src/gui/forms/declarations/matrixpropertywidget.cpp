#include "matrixpropertywidget.h"
#include "ui_matrixpropertywidget.h"

MatrixPropertyWidget::MatrixPropertyWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MatrixPropertyWidget)
{
    ui->setupUi(this);
}

MatrixPropertyWidget::~MatrixPropertyWidget()
{
    delete ui;
}
