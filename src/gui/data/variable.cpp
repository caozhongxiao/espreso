#include "variable.h"

Variable::Variable()
{
    this->mData = new ConstantType("");
}

Variable::Variable(const Variable& other) : NamedEntity(other)
{
    this->mData = other.mData->copy();
}

Variable::Variable(const QString& name, const DataType* data) : NamedEntity(name)
{
    this->mData = data->copy();
}

Variable::~Variable()
{
    delete this->mData;
}

QString Variable::toString() const
{
    return QString("%1 (%2)").arg(this->mName).arg(this->mData->toString());
}

DataType* Variable::data() const
{
    return this->mData;
}
