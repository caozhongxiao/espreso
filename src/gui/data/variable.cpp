#include "variable.h"

Variable::Variable()
{
    this->mName = "";
    this->mData = new StringType("");
}

Variable::Variable(const Variable &other)
{
    this->mData = other.mData->copy();
    this->mName = other.mName;
}

Variable::Variable(const QString& name, DataType* data)
{
    this->mName = name;
    this->mData = data;
}

Variable::~Variable()
{
}

QString Variable::toString() const
{
    return QString("%1 (%2)").arg(this->mName).arg(this->mData->toString());
}

QString Variable::name() const
{
    return this->mName;
}

DataType* Variable::data() const
{
    return this->mData;
}
