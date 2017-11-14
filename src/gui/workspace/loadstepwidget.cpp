#include "loadstepwidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

#include "../validators/validatorfactory.h"
#include "../elements/spinnerhandler.h"

using namespace espreso;

LoadstepWidget::LoadstepWidget(size_t id, Mesh* mesh, ECFObject* physics, QWidget* parent) :
    ECFObjectWidget(physics, parent)
{
    this->m_physics = physics;
    this->m_loadstep = m_physics->getParameter("load_steps_settings")->getParameter(QString::number(id).toStdString());
    this->m_obj = static_cast<ECFObject*>(m_loadstep);
    this->m_mesh = mesh;
    this->m_properties = nullptr;
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
    if (obj->name.compare("material_set") == 0
            || obj->name.compare("load_steps_settings") == 0)
        return;

    if ( obj->metadata.datatype.size() == 2 || obj->metadata.description.size() == 2)
    {
        if (this->m_properties == nullptr)
        {
            this->m_properties = new RegionPropertyWidget(m_mesh,
                                                          static_cast<PhysicsConfiguration*>(m_physics),
                                                          this->m_container,
                                                          tr("Boundary conditions"));
        }
        this->m_properties->addProperty(obj);
        this->m_widget->layout()->addWidget(m_properties);
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
