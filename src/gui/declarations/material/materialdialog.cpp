#include "materialdialog.h"
#include "ui_materialdialog.h"

#include <QComboBox>
#include <QMessageBox>
#include <QLabel>
#include <QScrollArea>
#include <QLineEdit>
#include <QDebug>
#include "../../elements/optionhandler.h"
#include "../../elements/formwidget.h"
#include "materialpropertytablewidget.h"

using namespace espreso;

MaterialDialog::MaterialDialog(MaterialConfiguration* material,
                               const QVector<std::string>& materialNames,
                               QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MaterialDialog)
{
    ui->setupUi(this);

    this->m_names = materialNames;
    this->m_material = material;
    this->drawMe();
}

MaterialDialog::~MaterialDialog()
{
    delete ui;
}

void MaterialDialog::drawMe()
{
    this->m_save.clear();
    this->m_valid.clear();

    this->m_frame = new QFrame(ui->frame);
    this->m_frameLayout = new QVBoxLayout;
    this->m_frame->setLayout(m_frameLayout);

    ui->frameLayout->addWidget(m_frame);

	this->iterateObject(m_material, this->m_frame);
}

void MaterialDialog::iterateObject(ECFObject* obj)
{
    QScrollArea* area = new QScrollArea(m_frame);
    QWidget* scrollWidget = new QWidget;
    QVBoxLayout* layout = new QVBoxLayout;
    scrollWidget->setLayout(layout);

    this->drawHeadline(obj, scrollWidget);

    this->processParameters(obj, scrollWidget);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);

    area->setWidgetResizable(true);
    area->setWidget(scrollWidget);
    this->m_frameLayout->addWidget(area);
}

void MaterialDialog::iterateObject(ECFObject* obj, QWidget* parent)
{
	// Parent widget determination
    QWidget* widget = parent;

    this->drawHeadline(obj, widget);

    this->processParameters(obj, widget);
}

void MaterialDialog::processParameters(ECFObject* obj, QWidget* widget)
{
    // Scalar properties
    MaterialPropertyTableWidget* propertyTable;
    bool propertyTableNotInserted = true;

    // Material details (name, description, ...)
    FormWidget* formWidget;
    bool formWidgetNotInserted = true;

    // Iterating over material parameters and creating proper UI widgets
    for (auto parameter = obj->parameters.cbegin();
         parameter != obj->parameters.cend();
         ++parameter)
    {
        if (!(*parameter)->metadata.isallowed())
            continue;

        if ( (*parameter)->isObject() )
        {
            this->iterateObject(static_cast<ECFObject*>(*parameter));
        }
        else if ((*parameter)->metadata.datatype.size())
        {

            ECFDataType type = (*parameter)->metadata.datatype.at(0);

            if ( type == ECFDataType::OPTION
                 || type == ECFDataType::ENUM_FLAGS )
            {
                this->drawOption(*parameter, widget);
            }
            else if ( type == ECFDataType::EXPRESSION )
            {

                if (propertyTableNotInserted)
                {
                    propertyTable = new MaterialPropertyTableWidget(widget);
                    this->m_save.append(propertyTable);
                    this->m_valid.append(propertyTable);
                    widget->layout()->addWidget(propertyTable);
                    propertyTableNotInserted = false;
                }
                propertyTable->addProperty(*parameter);

            }
            else if (type == ECFDataType::STRING)
            {
                if (formWidgetNotInserted)
                {
                    formWidget = new FormWidget(widget);
                    this->m_save.append(formWidget);
                    this->m_valid.append(formWidget);
                    widget->layout()->addWidget(formWidget);
                    formWidgetNotInserted = false;
                }
                formWidget->appendString(*parameter);
            }
        }
    }
}

void MaterialDialog::drawHeadline(ECFObject* obj, QWidget* widget)
{
    if (!obj->metadata.description.size())
        return;

    QString lblNameText = QString::fromStdString(obj->metadata.description.at(0));
    QLabel* lblName = new QLabel(lblNameText,
                                 widget);
    widget->layout()->addWidget(lblName);
}

void MaterialDialog::drawOption(ECFParameter* option, QWidget* widget)
{
    OptionHandler* handler = new OptionHandler(option, widget);
    widget->layout()->addWidget(handler);
    connect(handler, &OptionHandler::optionChanged, this, &MaterialDialog::redraw);
}

void MaterialDialog::redraw()
{
    foreach (ISavableObject* obj, m_save) {
        obj->save();
    }
    QFrame* tmp = this->m_frame;
    this->m_frame = nullptr;
    this->m_frameLayout = nullptr;
    ui->frameLayout->removeWidget(tmp);
    tmp->setParent(nullptr);
    this->drawMe();
}

void MaterialDialog::accept()
{
    foreach (IValidatableObject* obj, m_valid) {
        if (!obj->isValid())
        {
            QMessageBox msg;
            msg.setWindowTitle(tr("Error"));
            msg.setText(obj->errorMessage());
            msg.exec();
            return;
        }
    }


    ISavableObject* info = m_save.first();
    info->saveState();
    info->save();

    std::string newName =
            QString::fromStdString(
                m_material->getParameter(std::string("name"))->getValue()
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

        return;
    }

    foreach (ISavableObject* obj, m_save) {
        obj->save();
    }

    QDialog::accept();
}
