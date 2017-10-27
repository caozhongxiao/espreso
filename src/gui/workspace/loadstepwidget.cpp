#include "loadstepwidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

#include "../validators/validatorfactory.h"
#include "../elements/spinnerhandler.h"

using namespace espreso;

LoadstepWidget::LoadstepWidget(size_t id, ECFObject* physics, QWidget* parent) :
    ECFObjectWidget(physics, parent)
{
    this->m_physics = physics;
    this->m_loadstep = m_physics->getParameter("load_steps_settings")->getParameter(QString::number(id).toStdString());
    this->m_obj = static_cast<ECFObject*>(m_loadstep);
}


QWidget* LoadstepWidget::initContainer()
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

void LoadstepWidget::drawObject(ECFObject* obj)
{
    QWidget* widget = new QWidget(this->m_container);
    QLayout* layout;
    if (obj->parameters.size()) layout = new QFormLayout;
    else layout = new QVBoxLayout;
    widget->setLayout(layout);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);

    this->m_widget->layout()->addWidget(widget);


    if (obj->parameters.size())
    {
        this->processParameters(obj, widget);
    }
    else if ( !obj->parameters.size() && (obj->metadata.datatype.size() == 2) )
    {
        //TODO: INITIAL TEMPERATURE AND THICKNESS
    }
}

void LoadstepWidget::processParameters(ECFObject *obj, QWidget *widget)
{
    QFormLayout* l_layout = (QFormLayout*)widget->layout();

    for (auto parameter = obj->parameters.cbegin();
         parameter != obj->parameters.cend();
         ++parameter)
    {
        if (!(*parameter)->metadata.isallowed())
            continue;

        if ( (*parameter)->isObject() )
        {
            ECFObject* v_obj = static_cast<ECFObject*>(*parameter);
            if (v_obj->name.compare("material_set") == 0
                    || v_obj->name.compare("load_steps_settings") == 0)
                continue;
            this->drawObject(v_obj);
        }
        else if ((*parameter)->metadata.datatype.size() == 1)
        {
            ECFDataType type = (*parameter)->metadata.datatype.at(0);

            if ( type == ECFDataType::OPTION
                 || type == ECFDataType::ENUM_FLAGS )
            {
                l_layout->addRow(QString::fromStdString((*parameter)->metadata.description[0]),
                                 this->createOption(*parameter, widget, false));
            }
            else if ( type == ECFDataType::FLOAT )
            {
                DoubleValidatorFactory df;
                FieldHandler* handler = new FieldHandler(*parameter, &df);
                l_layout->addRow(QString::fromStdString((*parameter)->metadata.description[0]),
                        handler);
                this->m_savables.append(handler);
                this->m_validatables.append(handler);
            }
            else if ( type == ECFDataType::POSITIVE_INTEGER)
            {
                PositiveIntegerValidatorFactory df;
                FieldHandler* handler = new FieldHandler(*parameter, &df, false, widget);
                l_layout->addRow(QString::fromStdString((*parameter)->metadata.description[0]),
                        handler);
                this->m_savables.append(handler);
                this->m_validatables.append(handler);
            }
        }
    }
}
