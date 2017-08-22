#include "validatorfactory.h"

#include "expressionvalidator.h"
#include "../../data/common.h"


QValidator* ExpressionValidatorFactory::create(QObject *parent)
{
    return new ExpressionValidator(parent);
}


QValidator* NumberValidatorFactory::create(QObject *parent)
{
    return new QRegExpValidator(QRegExp(REGEXPR_DOUBLE), parent);
}
