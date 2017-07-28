#include "variable.h"

Variable::Variable()
{
    this->name = "";
    this->data = new StringType("");
}

Variable::Variable(const Variable &other)
{
    this->data = other.data->copy();
    this->name = other.name;
}

Variable::Variable(QString name, DataType* data)
{
    this->name = name;
    this->data = data;
}

Variable::~Variable()
{
}

QString Variable::toString() const
{
    return QString("%1 (%2)").arg(this->name).arg(this->data->toString());
}

QString Variable::name() const
{
    return this->name;
}

const DataType* Variable::data() const
{
    return this->data;
}
