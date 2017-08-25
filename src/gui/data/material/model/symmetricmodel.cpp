#include "symmetricmodel.h"

SymmetricModel::SymmetricModel() : PropertyModel()
{
    this->init();
}

SymmetricModel::SymmetricModel(const QList<DataType*>& values) :
    PropertyModel(values)
{
    if (values.size() != 6)
    {
        this->reportError();
        this->init();
        return;
    }

    foreach (DataType* val, values) {
        this->mMatrix << val;
    }
}

void SymmetricModel::setMatrix(const QList<DataType*>& values)
{
    if (values.size() != 6)
    {
        this->reportError();
        return;
    }

    this->clearMatrix();
    this->mMatrix = values;
}

void SymmetricModel::init()
{
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
}

