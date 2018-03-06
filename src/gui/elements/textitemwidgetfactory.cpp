#include "textitemwidgetfactory.h"

using namespace espreso;

TextItemWidgetFactory::TextItemWidgetFactory()
{

}


TextItemWidget* FilepathWidgetFactory::create(QWidget *parent)
{
    return new FilepathWidget(parent);
}


OptionWidgetFactory::OptionWidgetFactory(const QStringList &options)
{
    this->m_options = options;
}

TextItemWidget* OptionWidgetFactory::create(QWidget *parent)
{
    return new OptionWidget(this->m_options, parent);
}
