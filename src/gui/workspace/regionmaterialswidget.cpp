#include "regionmaterialswidget.h"
#include "ui_regionmaterialswidget.h"

using namespace espreso;

RegionMaterialsWidget::RegionMaterialsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::RegionMaterialsWidget)
{
    ui->setupUi(this);
}

RegionMaterialsWidget::~RegionMaterialsWidget()
{
    delete ui;
}
