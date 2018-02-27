#include "inputwidget.h"

#include "../elements/maptablewidgetfactory.h"
#include "../elements/textitemdelegate.h"
#include "../elements/textitemwidgetfactory.h"

using namespace espreso;

InputWidget::InputWidget(ECFObject* obj, QWidget* parent) :
    FixedECFObjectWidget(obj, parent)
{
}

void InputWidget::drawObject(ECFObject* obj)
{
    if (obj->metadata.datatype.size() == 2)
    {
        QWidget* widget = new QWidget;
        QLayout* layout = new QVBoxLayout;
        widget->setLayout(layout);

        this->createHeadline(obj, widget);

        MapTableWidgetFactory factory;
        MapTableWidget* mapwidget = factory.create(obj);
        this->m_savables.append(mapwidget);
        this->m_validatables.append(mapwidget);
        layout->addWidget(mapwidget);

        this->m_widget->layout()->addWidget(widget);

        return;
    }

    FixedECFObjectWidget::drawObject(obj);
}

ECFValueTableWidget* InputWidget::processString(ECFParameter* param, ECFValueTableWidget* table, QWidget* widget)
{
    ECFValueTableWidget* tw = this->createTableWidget(widget, table);
    if (param->name.compare("path") == 0)
    {
        tw->addWithDelegate(static_cast<ECFValue*>(param),
                            new TextItemDelegate(
                                QString::fromStdString(param->getValue()),
                                new FilepathWidgetFactory
                                ));
        return tw;
    }

    return FixedECFObjectWidget::processString(param, table, widget);
}
