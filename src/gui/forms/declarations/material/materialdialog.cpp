#include "materialdialog.h"
#include "ui_materialdialog.h"

#include <QComboBox>
#include "../../elements/optionhandler.h"
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

    this->iterateObject(m_material);
}

void MaterialDialog::iterateObject(ECFObject* obj)
{
    QFrame* frame = new QFrame(this->m_frame);
    frame->setFrameShape(QFrame::StyledPanel);
    QVBoxLayout* layout = new QVBoxLayout;
    frame->setLayout(layout);
    this->m_frameLayout->addWidget(frame);

    MaterialPropertyTableWidget* propertyTable;
    bool propertyTableNotInserted = true;

    for (auto parameter = obj->parameters.cbegin();
         parameter != obj->parameters.cend();
         ++parameter)
    {
        if (!(*parameter)->metadata.isallowed())
            continue;

        if ((*parameter)->metadata.datatype.size())
        {
            ECFDataType type = (*parameter)->metadata.datatype.at(0);
            if ( type == ECFDataType::OPTION )
            {
                this->drawOption(*parameter, frame);
            }
            else if ( type == ECFDataType::EXPRESSION )
            {
                if (propertyTableNotInserted)
                {
                    propertyTable = new MaterialPropertyTableWidget(frame);
                    layout->addWidget(propertyTable);
                    propertyTableNotInserted = false;
                }
                propertyTable->addProperty(*(*parameter));
            }
        }
        else if ( (*parameter)->isObject() )
        {
            this->iterateObject(static_cast<ECFObject*>(*parameter));
        }
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
    QFrame* tmp = this->m_frame;
    this->m_frame = nullptr;
    this->m_frameLayout = nullptr;
    ui->frameLayout->removeWidget(tmp);
    tmp->setParent(nullptr);
    this->drawMe();
}
