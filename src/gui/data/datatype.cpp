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
