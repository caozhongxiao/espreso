#include "regionobjectwidget.h"

using namespace espreso;

RegionObjectWidget::RegionObjectWidget(ECFObject* obj, QWidget *parent) :
    ECFObjectWidget(obj, parent)
{

}

QWidget* RegionObjectWidget::initContainer()
{
    QWidget* widget = new QWidget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    widget->setLayout(w_layout);

    return widget;
}

void RegionObjectWidget::drawObject(ECFObject* obj)
{
    QWidget* widget = new QWidget;
    QLayout* layout = new QVBoxLayout;
    widget->setLayout(layout);

    this->m_container->layout()->addWidget(widget);

    this->createHeadline(obj, widget);

    this->processParameters(obj, widget);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);
}
