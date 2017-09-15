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

MaterialDialog::MaterialDialog(ECFObject* material, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MaterialDialog)
{
    ui->setupUi(this);

    this->m_material = material;
    this->drawMe();
}

MaterialDialog::~MaterialDialog()
{
    delete ui;
}

void MaterialDialog::drawMe()
{
    this->m_frame = new QFrame(ui->frame);
    this->m_frameLayout = new QVBoxLayout;
    this->m_frame->setLayout(m_frameLayout);

    ui->frameLayout->addWidget(m_frame);

	this->iterateObject(m_material, this->m_frame);
}

void MaterialDialog::iterateObject(ECFObject* obj, QWidget* parent)
{
	// Parent widget determination
	QWidget* widget;
	QScrollArea* area = nullptr;
	if (parent == nullptr)
	{
		area = new QScrollArea(m_frame);
		QWidget* scrollWidget = new QWidget;
		QVBoxLayout* layout = new QVBoxLayout;
		scrollWidget->setLayout(layout);
		widget = scrollWidget;
	}
	else
	{
		widget = parent;
	}

	// Object name label
	if (obj->metadata.description.size())
	{
		QString lblNameText = QString::fromStdString(obj->metadata.description.at(0));
		QLabel* lblName = new QLabel(lblNameText,
									 widget);
		widget->layout()->addWidget(lblName);
	}

	// Scalar properties and object details (name, desc,...)
    MaterialPropertyTableWidget* propertyTable;
    bool propertyTableNotInserted = true;
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

	// ScrollArea should be drawn after its content is complete (according to Qt API)
	if (area != nullptr)
	{
		area->setWidgetResizable(true);
		area->setWidget(widget);
		this->m_frameLayout->addWidget(area);
	}
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

    foreach (ISavableObject* obj, m_save) {
        obj->save();
    }

    QDialog::accept();
}
