#ifndef DATATYPE_H
#define DATATYPE_H

#include <QString>

class DataType
{

public:
    DataType();
    virtual ~DataType() {}
    virtual QString toString() const = 0;
    virtual DataType* copy() const = 0;
};

class StringType : public DataType
{
private:
    QString data;

public:
    StringType(QString data);
    QString toString() const override;
    DataType* copy() const override;
};

#endif // DATATYPE_H
