#include "materialwidget.h"

#include <QVBoxLayout>
#include <QScrollArea>
#include <QComboBox>
#include <QMessageBox>
#include <QLabel>
#include <QScrollArea>
#include <QLineEdit>
#include <QDebug>
#include "../../elements/optionhandler.h"
#include "../../elements/formwidget.h"
#include "../../elements/boolhandler.h"
#include "materialpropertytablewidget.h"

using namespace espreso;

MaterialWidget::MaterialWidget(MaterialConfiguration* material,
                               const QVector<std::string>& materialNames,
                               QWidget *parent) :
    ECFObjectWidget(material, parent)
{
    this->m_names = materialNames;
}


QWidget* MaterialWidget::initContainer()
{
    QScrollArea* area = new QScrollArea;

    QWidget* widget = new QWidget(area);
    this->m_widget = widget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    widget->setLayout(w_layout);

    this->drawHeadline(this->m_obj, widget);

    area->setWidgetResizable(true);
    area->setWidget(widget);

    return area;
}

void MaterialWidget::drawObject(ECFObject* obj)
{
    QWidget* widget = new QWidget;
    QVBoxLayout* layout = new QVBoxLayout;
    widget->setLayout(layout);

    this->m_widget->layout()->addWidget(widget);

    this->drawHeadline(obj, widget);

    this->processParameters(obj, widget);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);
}

void MaterialWidget::drawHeadline(ECFObject* obj, QWidget* widget)
{
    if (!obj->metadata.description.size())
        return;

    QString lblNameText = QString::fromStdString(obj->metadata.description.at(0));
    QLabel* lblName = new QLabel(lblNameText);
    widget->layout()->addWidget(lblName);
}

void MaterialWidget::processParameters(ECFObject *obj, QWidget *widget)
{
    // Scalar properties
    MaterialPropertyTableWidget* propertyTable;
    bool propertyTableNotInserted = true;

    // Material details (name, description, ...)
    FormWidget* formWidget = nullptr;

    // Iterating over material parameters and creating proper UI widgets
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
        else if ((*parameter)->metadata.datatype.size())
        {

            ECFDataType type = (*parameter)->metadata.datatype.at(0);

            if ( type == ECFDataType::OPTION
                 || type == ECFDataType::ENUM_FLAGS )
            {
                QWidget* handler = this->createOption(*parameter, widget);
                widget->layout()->addWidget(handler);
            }
            else if (type == ECFDataType::BOOL)
            {
                QWidget* handler = this->createBool(*parameter, widget);
                widget->layout()->addWidget(handler);
            }
            else if ( type == ECFDataType::EXPRESSION )
            {

                if (propertyTableNotInserted)
                {
                    propertyTable = new MaterialPropertyTableWidget(widget);
                    this->m_savables.append(propertyTable);
                    this->m_validatables.append(propertyTable);
                    widget->layout()->addWidget(propertyTable);
                    propertyTableNotInserted = false;
                }
                propertyTable->addProperty(*parameter);

            }
            else if (type == ECFDataType::STRING)
            {
                formWidget = this->createFormWidget(widget, formWidget);
                formWidget->appendString(*parameter);
            }
            else if (type == ECFDataType::FLOAT)
            {
                formWidget = this->createFormWidget(widget, formWidget);
                formWidget->appendFloat(*parameter);
            }
            else if (type == ECFDataType::NONNEGATIVE_INTEGER)
            {
                formWidget = this->createFormWidget(widget, formWidget);
                formWidget->appendNonnegativeInteger(*parameter);
            }
            else if (type == ECFDataType::POSITIVE_INTEGER)
            {
                formWidget = this->createFormWidget(widget, formWidget);
                formWidget->appendPositiveInteger(*parameter);
            }
        }
    }
}

FormWidget* MaterialWidget::createFormWidget(QWidget* container, FormWidget* form)
{
    if (form != nullptr) return form;

    FormWidget* formWidget = new FormWidget;
    this->m_savables.append(formWidget);
    this->m_validatables.append(formWidget);
    container->layout()->addWidget(formWidget);

    return formWidget;
}

bool MaterialWidget::isValid()
{
    if (!ECFObjectWidget::isValid())
        return false;

    ISavableObject* info = m_savables.first();
    info->saveState();
    info->save();

    std::string newName =
            QString::fromStdString(
                m_obj->getParameter(std::string("name"))->getValue()
                )
                .toUpper()
                .toStdString();

    if (m_names.indexOf(newName) != -1)
    {
        QMessageBox msg;
        msg.setWindowTitle(tr("Error"));
        msg.setText(tr("Material with same name already exists!"));
        msg.exec();

        info->restoreState();

        return false;
    }

    foreach (ISavableObject* obj, m_savables) {
        obj->save();
    }

    return true;
}
