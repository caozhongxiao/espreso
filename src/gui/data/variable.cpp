#include "variable.h"

Variable::Variable()
{
    this->mData = new ConstantType("");
}

//Variable::Variable(const Variable& other)
//{
//    this->mData = other.mData->copy();
//    this->mName = other.mName;
//}

Variable::Variable(const QString& name, DataType* data) : NamedEntity(name)
{
    this->mData = data;
}

//Variable::~Variable()
//{
//}

QString Variable::toString() const
{
    return QString("%1 (%2)").arg(this->mName).arg(this->mData->toString());
}

DataType* Variable::data() const
{
    return this->mData;
}
