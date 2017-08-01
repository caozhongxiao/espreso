#ifndef DATATYPE_H
#define DATATYPE_H

#include <QString>
#include <QVariant>

enum DTLib
{
    CONSTANT                = 0,
    FUNCTION                = 1,
    TABLE                   = 2,
    PIECEWISE_FUNCTION      = 3,
    STRING                  = 4
};

class DataType
{

public:
    DataType();
    virtual ~DataType() {}
    virtual QString toString() const = 0;
    virtual DataType* copy() const = 0;
    virtual int type() const = 0;
};

class StringType : public DataType
{
protected:
    QString data;

public:
    StringType(const QString& data);
    QString toString() const override;
    DataType* copy() const override;
    virtual int type() const override;
};

class ConstantType : public StringType
{
public:
    ConstantType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    int type() const override;
};

class FunctionType : public StringType
{
public:
    FunctionType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    int type() const override;
};

#endif // DATATYPE_H
