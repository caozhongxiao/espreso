#include "validatorfactory.h"

#include "../data/common.h"

QValidator* DoubleValidatorFactory::create(QObject *parent) const
{
    return new QRegExpValidator(QRegExp(GUI_REGEXPR_DOUBLE), parent);
}

QValidator* PositiveIntegerValidatorFactory::create(QObject *parent) const
{
    return new QRegExpValidator(QRegExp("^[0]?([1-9][0-9]*)$"), parent);
}
