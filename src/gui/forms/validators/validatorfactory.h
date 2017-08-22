#ifndef VALIDATORFACTORY_H
#define VALIDATORFACTORY_H

#include <QValidator>

class ValidatorFactory
{

public:
    ValidatorFactory() {}
    virtual ~ValidatorFactory() {};

    virtual QValidator* create(QObject* parent = nullptr) = 0;
};

class ExpressionValidatorFactory : public ValidatorFactory
{
public:
    QValidator* create(QObject *parent = nullptr) override;
};

class NumberValidatorFactory : public ValidatorFactory
{
public:
    QValidator* create(QObject *parent = nullptr) override;
};

#endif // VALIDATORFACTORY_H
