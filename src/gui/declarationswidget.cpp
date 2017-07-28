#include "declarationswidget.h"
#include "ui_declarationswidget.h"

DeclarationsWidget::DeclarationsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DeclarationsWidget)
{
    ui->setupUi(this);
}

DeclarationsWidget::~DeclarationsWidget()
{
    delete ui;
}
