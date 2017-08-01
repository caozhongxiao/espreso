#include "datatype.h"

DataType::DataType()
{

}

StringType::StringType(const QString& data):DataType()
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


int ConstantType::type() const
{
    return DTLib::CONSTANT;
}

DataType* ConstantType::copy() const
{
    return new ConstantType(this->data);
}


int FunctionType::type() const
{
    return DTLib::FUNCTION;
}

DataType* FunctionType::copy() const
{
    return new FunctionType(this->data);
}
