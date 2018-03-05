#include "textitemwidgetfactory.h"

using namespace espreso;

TextItemWidgetFactory::TextItemWidgetFactory()
{

}


TextItemWidget* FilepathWidgetFactory::create(QWidget *parent)
{
    return new FilepathWidget(parent);
}
