#include "outputconfigurationwidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

#include "../elements/ecfobjectwidgetfactory.h"
#include "../elements/integertabwidget.h"
#include <memory>

using namespace espreso;

OutputConfigurationWidget::OutputConfigurationWidget(ECFObject* output, QWidget* parent) :
    ECFObjectWidget(output, parent)
{

}

QWidget* OutputConfigurationWidget::initContainer()
{
    QScrollArea* area = new QScrollArea;

    QWidget* widget = new QWidget(area);
    this->m_widget = widget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    widget->setLayout(w_layout);

    area->setWidgetResizable(true);
    area->setWidget(widget);

    return area;
}

void OutputConfigurationWidget::drawObject(ECFObject* obj)
{
    if (obj->name.compare("monitoring") == 0)
    {
        std::unique_ptr<ECFObjectWidgetFactory> factory(new FixedECFObjectWidgetFactory(false));
        this->m_widget->layout()->addWidget(new IntegerTabWidget(obj, std::move(factory)));
        return;
    }

    QWidget* widget = new QWidget(this->m_container);
    QLayout* layout = new QVBoxLayout;
    widget->setLayout(layout);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);

    this->m_widget->layout()->addWidget(widget);

    this->createHeadline(obj, widget);

    this->processParameters(obj, widget);
}
