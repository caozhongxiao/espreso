#include "optionhandler.h"

using namespace espreso;

#include <QVBoxLayout>
#include <QComboBox>

OptionHandler::OptionHandler(ECFParameter* option, QWidget* parent) :
    QWidget(parent)
{
    this->m_option = option;

    QVBoxLayout* layout = new QVBoxLayout;
    this->setLayout(layout);

    QComboBox* cmb = new QComboBox(this);
    layout->addWidget(cmb);

    int index = 0;
    for (auto item = option->metadata.options.cbegin();
         item != option->metadata.options.cend()
         && item->isallowed();
         ++item)
    {
        cmb->addItem(QString::fromStdString( item->name ));
        if (item->name.compare(option->getValue()) == 0)
        {
            cmb->setCurrentIndex(index);
        }
        index++;
    }

    this->optionsAdded = true;

    connect(cmb, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &OptionHandler::onIndexChanged);
}

void OptionHandler::onIndexChanged(int index)
{
    if (!optionsAdded)
        return;

    this->m_option->setValue(m_option->metadata.options.at(index).name);

    emit optionChanged();
}
