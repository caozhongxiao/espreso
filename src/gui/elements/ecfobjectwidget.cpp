#include "ecfobjectwidget.h"
#include "ui_ecfobjectwidget.h"

#include <QMessageBox>
#include <QScrollArea>
#include <QLabel>

using namespace espreso;

ECFObjectWidget::ECFObjectWidget(ECFObject* object, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ECFObjectWidget)
{
    this->m_obj = object;
    ui->setupUi(this);
}

ECFObjectWidget::~ECFObjectWidget()
{
    delete ui;
}

void ECFObjectWidget::init()
{
    this->drawMe();
}

void ECFObjectWidget::drawMe()
{
    this->m_savables.clear();
    this->m_validatables.clear();

    this->m_container = this->initContainer();
    ui->layout->addWidget(this->m_container);

    this->drawObject(this->m_obj);
}

OptionHandler* ECFObjectWidget::createOption(ECFParameter* option, QWidget* parent, bool withLabel)
{
    OptionHandler* handler = new OptionHandler(option, parent, withLabel);
    connect(handler, &OptionHandler::optionChanged, this, &ECFObjectWidget::redraw);

    return handler;
}

BoolHandler* ECFObjectWidget::createBool(ECFParameter* param, QWidget* parent)
{
    BoolHandler* handler = new BoolHandler(param, parent);
    connect(handler, &BoolHandler::stateChanged, this, &ECFObjectWidget::redraw);

    return handler;
}

void ECFObjectWidget::createHeadline(ECFObject* obj, QWidget* widget)
{
    if (!obj->metadata.description.size())
        return;

    QString lblNameText = QString::fromStdString(obj->metadata.description.at(0));
    QLabel* lblName = new QLabel(lblNameText);
    widget->layout()->addWidget(lblName);
}

void ECFObjectWidget::redraw()
{
    foreach (ISavableObject* obj, m_savables) {
        obj->save();
    }

    ui->layout->removeWidget(this->m_container);

    this->drawMe();
}

bool ECFObjectWidget::validate()
{
    foreach (IValidatableObject* obj, m_validatables) {
        if (!obj->isValid())
        {
            QMessageBox msg;
            msg.setWindowTitle(tr("Error"));
            msg.setText(obj->errorMessage());
            msg.exec();
            return false;
        }
    }

    return true;
}

bool ECFObjectWidget::isValid()
{
    foreach (IValidatableObject* obj, m_validatables) {
        if (!obj->isValid())
        {
            this->m_errormsg = obj->errorMessage();
            return false;
        }
    }

    return true;
}

QString ECFObjectWidget::errorMessage()
{
    this->isValid();

    return this->m_errormsg;
}

void ECFObjectWidget::save()
{
    foreach (ISavableObject* obj, m_savables)
    {
        obj->save();
    }
}

void ECFObjectWidget::processParameters(ECFObject* obj, QWidget* widget)
{
    // Property table (name, expression, units,...)
    MaterialPropertyTableWidget* propertyTable = nullptr;
    // Widget with form layout. Container for form fields.
    FormWidget* formWidget = nullptr;

    // Iterating object parameters and creating UI elements for them
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
                formWidget = this->processOptionEnum(*parameter, formWidget, widget);
            }
            else if (type == ECFDataType::BOOL)
            {
                formWidget = this->processBool(*parameter, formWidget, widget);
            }
            else if (type == ECFDataType::EXPRESSION)
            {
                propertyTable = this->processExpression(*parameter, propertyTable, widget);
            }
            else if (type == ECFDataType::STRING)
            {
                formWidget = this->processString(*parameter, formWidget, widget);
            }
            else if (type == ECFDataType::FLOAT)
            {
                formWidget = this->processFloat(*parameter, formWidget, widget);
            }
            else if (type == ECFDataType::NONNEGATIVE_INTEGER)
            {
                formWidget = this->processNonnegativeInteger(*parameter, formWidget, widget);
            }
            else if (type == ECFDataType::POSITIVE_INTEGER)
            {
                formWidget = this->processPositiveInteger(*parameter, formWidget, widget);
            }
            else if (type == ECFDataType::REGION)
            {
                formWidget = this->processRegion(*parameter, formWidget, widget);
            }
        }
    }
}

FormWidget* ECFObjectWidget::processOptionEnum(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    QWidget* handler = this->createOption(parameter, widget, false);
    formWidget->appendRow(QString::fromStdString(parameter->metadata.description[0]),
            handler);

    return formWidget;
}

FormWidget* ECFObjectWidget::processBool(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    QWidget* handler = this->createBool(parameter, widget);
    formWidget->appendRow(QString::fromStdString(parameter->metadata.description[0]),
            handler);

    return formWidget;
}

MaterialPropertyTableWidget* ECFObjectWidget::processExpression(ECFParameter* parameter, MaterialPropertyTableWidget* tableWidget, QWidget* widget)
{
    MaterialPropertyTableWidget* propertyTable = tableWidget;
    if (propertyTable == nullptr)
    {
        propertyTable = new MaterialPropertyTableWidget(widget);
        this->m_savables.append(propertyTable);
        this->m_validatables.append(propertyTable);
        widget->layout()->addWidget(propertyTable);
    }
    propertyTable->addProperty(parameter);

    return propertyTable;
}

FormWidget* ECFObjectWidget::processString(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    formWidget->appendString(parameter);

    return formWidget;
}

FormWidget* ECFObjectWidget::processFloat(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    formWidget->appendFloat(parameter);

    return formWidget;
}

FormWidget* ECFObjectWidget::processNonnegativeInteger(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    formWidget->appendNonnegativeInteger(parameter);

    return formWidget;
}

FormWidget* ECFObjectWidget::processPositiveInteger(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    formWidget->appendPositiveInteger(parameter);

    return formWidget;
}

FormWidget* ECFObjectWidget::processRegion(ECFParameter* parameter, FormWidget* form, QWidget* widget)
{
    FormWidget* formWidget = this->createFormWidget(widget, form);
    formWidget->appendString(parameter);

    return formWidget;
}

FormWidget* ECFObjectWidget::createFormWidget(QWidget* container, FormWidget* form)
{
    if (form != nullptr) return form;

    FormWidget* formWidget = new FormWidget;
    this->m_savables.append(formWidget);
    this->m_validatables.append(formWidget);
    container->layout()->addWidget(formWidget);

    return formWidget;
}
