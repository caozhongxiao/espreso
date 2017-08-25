#include "anisotropicmodel.h"

AnisotropicModel::AnisotropicModel() : PropertyModel()
{
    this->init();
}

AnisotropicModel::AnisotropicModel(const QList<DataType*>& values) :
    PropertyModel(values)
{
    if (values.size() != 9)
    {
        this->reportError();
        this->init();
        return;
    }

    foreach (DataType* val, values) {
        this->mMatrix << val;
    }
}

void AnisotropicModel::setMatrix(const QList<DataType*>& values)
{
    if (values.size() != 9)
    {
        this->reportError();
        return;
    }

    this->clearMatrix();
    this->mMatrix = values;
}

void AnisotropicModel::init()
{
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
}
