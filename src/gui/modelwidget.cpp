#include "modelwidget.h"
#include "ui_modelwidget.h"

ModelWidget::ModelWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ModelWidget)
{
    ui->setupUi(this);
}

ModelWidget::~ModelWidget()
{
    delete ui;
}
