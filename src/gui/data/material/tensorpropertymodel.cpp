#include "tensorpropertymodel.h"

TensorPropertyModel::TensorPropertyModel() :
    NamedEntity()
{
    this->mSize = 1;
    this->mMatrix.append(new DummyType());
}

TensorPropertyModel::TensorPropertyModel(int size, const QString& name) :
    NamedEntity(name)
{
    this->mSize = size;
    for (int i = 0; i < size * size; ++i)
    {
        this->mMatrix.append(new DummyType());
    }
}

TensorPropertyModel::TensorPropertyModel(int size, const QString& name, const QList<DataType*>& values)
    : TensorPropertyModel(size, name)
{
    if (values.size() != size*size)
    {
        qWarning("TensorPropertyModel: Size does not match the number of values in matrix!");
        return;
    }

    foreach (DataType* val, values) {
        this->mMatrix << val;
    }
}

TensorPropertyModel::TensorPropertyModel(const TensorPropertyModel& other) :
    NamedEntity(other)
{
    this->mSize = other.mSize;
    QVector<DataType*> tmp;
    foreach (DataType* data, other.mMatrix) {
        tmp << data->copy();
    }
    this->mMatrix = tmp;
}

TensorPropertyModel::~TensorPropertyModel()
{
    foreach (DataType* val, this->mMatrix) {
        delete val;
    }
    this->mMatrix.clear();
}

QVector<DataType*> TensorPropertyModel::matrix() const
{
    return this->mMatrix;
}

void TensorPropertyModel::clearMatrix()
{
    foreach (DataType* val, this->mMatrix) {
        delete val;
    }
    this->mMatrix.clear();
}

void TensorPropertyModel::setMatrix(const QList<DataType*>& values)
{
    if (values.size() != mSize * mSize)
    {
        qWarning("TensorPropertyModel: List size does not match the property size!");
        return;
    }

    this->clearMatrix();

    foreach (DataType* val, values) {
        this->mMatrix << val;
    }
}

void TensorPropertyModel::setCellValue(int row, int column, DataType* value)
{
    int index = row * mSize + column;
    if (index >= mSize * mSize)
    {
        qWarning("TensorPropertyModel: Index [%d][%d] out of range!", row, column);
        return;
    }

    DataType* tmp = this->mMatrix.at(index);
    this->mMatrix[index] = value;
    delete tmp;
}

DataType* TensorPropertyModel::cellValue(int row, int column) const
{
    int index = row * mSize + column;
    if (index >= mSize * mSize)
    {
        qWarning("TensorPropertyModel: Index [%d][%d] out of range!", row, column);
        return nullptr;
    }

    return this->mMatrix.at(index);
}

int TensorPropertyModel::size() const
{
    return this->mSize;
}
