#include "datatype.h"

DataType::DataType()
{

}

StringType::StringType(QString data):DataType()
{
    this->data = data;
}

QString StringType::toString() const
{
    return this->data;
}

DataType* StringType::copy() const
{
    return new StringType(this->data);
}

int StringType::type() const
{
    return DTLib::STRING;
}


ConstantType::ConstantType(const QString& data):DataType()
{
    this->data = data;
}

QString ConstantType::toString() const
{
    return this->data;
}

DataType* ConstantType::copy() const
{
    return new ConstantType(this->data);
}

int ConstantType::type() const
{
    return DTLib::CONSTANT;
}
