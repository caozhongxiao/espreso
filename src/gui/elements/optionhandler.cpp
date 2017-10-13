#include "optionhandler.h"
#include "ui_optionhandler.h"

using namespace espreso;


OptionHandler::OptionHandler(ECFParameter* option, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::OptionHandler)
{
    ui->setupUi(this);

    this->m_option = option;

    if (option->metadata.description.size())
    {
        ui->label->setText(QString::fromStdString(option->metadata.description.at(0)));
    }

    int index = 0;
    for (auto item = option->metadata.options.cbegin();
         item != option->metadata.options.cend()
         && item->isallowed();
         ++item)
    {
        ui->cmb->addItem(QString::fromStdString( item->name ));
        if (item->name.compare(option->getValue()) == 0)
        {
            ui->cmb->setCurrentIndex(index);
        }
        index++;
    }

    this->optionsAdded = true;

    connect(ui->cmb, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &OptionHandler::onIndexChanged);
}

void OptionHandler::onIndexChanged(int index)
{
    if (!optionsAdded)
        return;

    this->m_option->setValue(m_option->metadata.options.at(index).name);

    emit optionChanged();
}

OptionHandler::~OptionHandler()
{
    delete ui;
}
