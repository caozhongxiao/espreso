#include "expressionvalidator.h"

#include "../../data/common.h"
#include "../../expression.h"

ExpressionValidator::ExpressionValidator(QObject* parent) :
    QValidator(parent)
{

}

QValidator::State ExpressionValidator::validate(QString& input, int& pos) const
{
    bool valid = Expression::isValid(input.toStdString(), Common::fnVariables());

    if (valid)
        return QValidator::Acceptable;

    valid = Expression::isValid(input.toStdString(), Common::fnVariablesParts());

    if (valid)
        return QValidator::Intermediate;

    return QValidator::Invalid;
}
