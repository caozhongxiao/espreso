#include "diagonalmodel.h"

DiagonalModel::DiagonalModel() : PropertyModel()
{
    this->init();
}

DiagonalModel::DiagonalModel(const QList<DataType*>& values) :
    PropertyModel(values)
{
    if (values.size() != 3)
    {
        this->reportError();
        this->init();
        return;
    }

    foreach (DataType* val, values) {
        this->mMatrix << val;
    }
}

void DiagonalModel::setMatrix(const QList<DataType*>& values)
{
    if (values.size() != 3)
    {
        this->reportError();
        return;
    }

    this->clearMatrix();
    this->mMatrix = values;
}

void DiagonalModel::init()
{
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
    this->mMatrix << new ExpressionType("0");
}
