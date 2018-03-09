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

TextWidgetFactory::TextWidgetFactory(ValidatorFactory *validatorFactory)
{
    this->m_factory = validatorFactory;
}

TextWidgetFactory::~TextWidgetFactory()
{
    if (this->m_factory != nullptr) delete this->m_factory;
}

TextItemWidget* TextWidgetFactory::create(QWidget *parent)
{
    TextWidget* widget = new TextWidget(parent);
    if (this->m_factory != nullptr)
        widget->setValidator(this->m_factory->create());
    return widget;
}
