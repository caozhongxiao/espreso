#include "physicswidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>
#include <QComboBox>

using namespace espreso;

PhysicsWidget::PhysicsWidget(ECFConfiguration* ecf, QWidget* parent) :
    ECFObjectWidget(ecf, parent)
{
    this->m_ecf = ecf;
}

QWidget* PhysicsWidget::initContainer()
{
    ECFParameter* physics = m_ecf->getParameter("physics");

    QScrollArea* area = new QScrollArea;

    QWidget* widget = new QWidget(area);
    this->m_widget = widget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    widget->setLayout(w_layout);

    QComboBox* cmbPhysics = new QComboBox(widget);
    w_layout->addWidget(cmbPhysics);

    int active = 0;
    int index = 0;
    for (auto option = physics->metadata.options.begin(); option != physics->metadata.options.end(); ++option, ++index)
    {
        QString name = QString::fromStdString(option->name);
        cmbPhysics->addItem(name);

        if (option->name.compare(physics->getValue()) == 0) active = index;
    }

    cmbPhysics->setCurrentIndex(active);
    this->m_obj = this->physics(active);
    connect(cmbPhysics, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &PhysicsWidget::onPhysicsChange);

    area->setWidgetResizable(true);
    area->setWidget(widget);

    return area;
}

ECFObject* PhysicsWidget::physics(int index)
{
    switch (index)
    {
        case 0:
            return &m_ecf->heat_transfer_2d;
        case 1:
            return &m_ecf->heat_transfer_3d;
        case 2:
            return &m_ecf->structural_mechanics_2d;
        case 3:
            return &m_ecf->structural_mechanics_3d;
        default:
            qCritical() << tr("WorkflowWidget: Invalid index of physics!");
            return nullptr;
    }
}

void PhysicsWidget::onPhysicsChange(int index)
{
    ECFParameter* physics = m_ecf->getParameter("physics");
    physics->setValue(physics->metadata.options[index].name);

    this->m_obj = this->physics(index);

    this->redraw();
}

void PhysicsWidget::drawObject(ECFObject* obj)
{
    QScrollArea* area = new QScrollArea(this->m_container);
    QWidget* scrollWidget = new QWidget;
    QFormLayout* layout = new QFormLayout;
    scrollWidget->setLayout(layout);

    this->processParameters(obj, scrollWidget);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);

    area->setWidgetResizable(true);
    area->setWidget(scrollWidget);
    this->m_widget->layout()->addWidget(area);
}

void PhysicsWidget::processParameters(ECFObject* obj, QWidget* widget)
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
            this->drawObject(static_cast<ECFObject*>(*parameter));
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

            }
            else if ( type == ECFDataType::POSITIVE_INTEGER)
            {

            }
            else
            {
                if ((*parameter)->metadata.description.size() > 0)
                    qInfo() << QString::fromStdString((*parameter)->metadata.description.at(0));
            }
        }
    }
}
