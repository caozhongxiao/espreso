#include "isotropicmodel.h"

IsotropicModel::IsotropicModel() : PropertyModel()
{
    this->init();
}

IsotropicModel::IsotropicModel(const QList<DataType*>& values) :
    PropertyModel(values)
{
    if (values.size() != 1)
    {
        this->reportError();
        this->init();
        return;
    }

    foreach (DataType* val, values) {
        this->mMatrix << val;
    }
}

void IsotropicModel::setMatrix(const QList<DataType*>& values)
{
    if (values.size() != 1)
    {
        this->reportError();
        return;
    }

    this->clearMatrix();
    this->mMatrix = values;
}


void IsotropicModel::init()
{
    this->mMatrix << new ExpressionType("0");
}
