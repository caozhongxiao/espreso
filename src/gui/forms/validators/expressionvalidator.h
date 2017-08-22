#ifndef EXPRESSIONVALIDATOR_H
#define EXPRESSIONVALIDATOR_H

#include <QValidator>

class ExpressionValidator : public QValidator
{
public:
    ExpressionValidator(QObject* parent = nullptr);
    QValidator::State validate(QString &, int &) const;
};

#endif // EXPRESSIONVALIDATOR_H
