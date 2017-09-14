#include "tensorpropertymodel.h"

TensorPropertyModel::TensorPropertyModel() :
    NamedEntity()
{
    this->mSize = 1;
    this->mMatrix.append(TensorPropertyModelItem());
}

TensorPropertyModel::TensorPropertyModel(int size, const QString& name,
                                         const QList<TensorPropertyModelItem>& values) :
    NamedEntity(name)
{
//    if (values.size() != size*size)
//    {
//        qWarning("TensorPropertyModel: Size does not match the number of values in matrix!");
//        return;
//    }

    foreach (TensorPropertyModelItem val, values) {
        this->mMatrix << val;
    }
}

TensorPropertyModel::TensorPropertyModel(const TensorPropertyModel& other) :
    NamedEntity(other)
{
    this->mSize = other.mSize;
    this->mMatrix = other.mMatrix;
}

TensorPropertyModel::~TensorPropertyModel()
{
    this->mMatrix.clear();
}

const QVector<TensorPropertyModelItem>& TensorPropertyModel::items() const
{
    return this->mMatrix;
}

//void TensorPropertyModel::setMatrix(const QList<TensorPropertyModelItem>& values)
//{
//    if (values.size() != mSize * mSize)
//    {
//        qWarning("TensorPropertyModel: List size does not match the property size!");
//        return;
//    }

//    this->mMatrix.clear();

//    foreach (TensorPropertyModelItem val, values) {
//        this->mMatrix << val;
//    }
//}

void TensorPropertyModel::setCellValue(int row, int column, TensorPropertyModelItem value)
{
//    int index = row * mSize + column;
//    if (index >= mSize * mSize)
//    {
//        qWarning("TensorPropertyModel: Index [%d][%d] out of range!", row, column);
//        return;
//    }

//    this->mMatrix[index] = value;
}

const TensorPropertyModelItem& TensorPropertyModel::cellValue(int row, int column) const
{
//    int index = row * mSize + column;
//    if (index >= mSize * mSize)
//    {
//        qCritical("TensorPropertyModel: Index [%d][%d] out of range!", row, column);
//    }

//    return this->mMatrix.at(index);
}

int TensorPropertyModel::size() const
{
    return this->mSize;
}
