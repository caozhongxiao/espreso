#include "inputwidget.h"

#include "../elements/maptablewidgetfactory.h"
#include "../elements/pathhandler.h"

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

FormWidget* InputWidget::processString(ECFParameter* param, FormWidget* fw, QWidget* widget)
{
    FormWidget* form = this->createFormWidget(widget, fw);
    if (param->name.compare("path") == 0)
    {
        PathHandler* ph = new PathHandler(param, this);
        this->m_validatables.append(ph);
        form->appendRow(
                QString::fromStdString(param->metadata.description[0]),
                ph);
        return form;
    }

    return FixedECFObjectWidget::processString(param, form, widget);
}
