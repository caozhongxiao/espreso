#include "propertymodel.h"

PropertyModel::PropertyModel()
{
}

PropertyModel::PropertyModel(const QList<DataType*>& values)
{
}

PropertyModel::~PropertyModel()
{
    foreach (DataType* val, this->mMatrix) {
        delete val;
    }
    this->mMatrix.clear();
}

QList<DataType*> PropertyModel::matrix() const
{
    return this->mMatrix;
}

void PropertyModel::reportError()
{
    qWarning("Invalid number of matrix values! New matrix was not used!");
}

void PropertyModel::clearMatrix()
{
    foreach (DataType* val, this->mMatrix) {
        delete val;
    }
    this->mMatrix.clear();
}
