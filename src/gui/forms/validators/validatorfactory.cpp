#include "validatorfactory.h"

#include "../../data/common.h"

QValidator* NumberValidatorFactory::create(QObject *parent)
{
    return new QRegExpValidator(QRegExp(GUI_REGEXPR_DOUBLE), parent);
}
