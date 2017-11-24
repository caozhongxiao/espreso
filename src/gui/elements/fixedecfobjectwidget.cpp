#include "fixedecfobjectwidget.h"

using namespace espreso;

FixedECFObjectWidget::FixedECFObjectWidget(ECFObject* obj, QWidget* parent) :
    ECFObjectWidget(obj, parent)
{
    this->m_draw_headlines = true;
}

QWidget* FixedECFObjectWidget::initContainer()
{
    QWidget* widget = new QWidget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    widget->setLayout(w_layout);

    return widget;
}

void FixedECFObjectWidget::drawObject(ECFObject* obj)
{
    QWidget* widget = new QWidget;
    QLayout* layout = new QVBoxLayout;
    widget->setLayout(layout);

    this->m_container->layout()->addWidget(widget);

    if (this->m_draw_headlines) this->createHeadline(obj, widget);

    this->processParameters(obj, widget);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);
}

void FixedECFObjectWidget::setDrawHeadline(bool draw)
{
    this->m_draw_headlines = draw;
}
