#ifndef DATATYPE_H
#define DATATYPE_H

#include <QString>
#include <QVariant>
#include <QVector>
#include <QPair>
#include <QObject>

//enum DTLib
//{
//    CONSTANT                = 0,
//    FUNCTION                = 1,
//    TABLE                   = 2,
//    PIECEWISE_FUNCTION      = 3,
//    STRING                  = 4,
//    VARIABLE_LINK           = 5
//};

class DataTypeVisitor;

class DataType
{

public:
    DataType() {}
    virtual ~DataType() {}
    virtual QString toString() const = 0;
    virtual DataType* copy() const = 0;
    virtual void accept(DataTypeVisitor* visitor) const = 0;
};

class StringType : public DataType
{
protected:
    QString data;
    StringType(const QString& data);

public:
    QString toString() const override;
    virtual DataType* copy() const = 0;
    virtual void accept(DataTypeVisitor* visitor) const = 0;
};

class ConstantType : public StringType
{
public:
    ConstantType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) const;
};

class FunctionType : public StringType
{
public:
    FunctionType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) const;
};

class VariableLinkType : public StringType
{
public:
    VariableLinkType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) const;
};

class TableType : public DataType
{
protected:
    QVector<QPair<QString, QString> > mRows;

public:
    TableType(const QVector<QPair<QString, QString> >& data);
    QString toString() const override;
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) const;

    QVector<QPair<QString, QString> > data() const;
};

class PiecewiseFunctionType : public DataType
{
protected:
    QVector<QVector<QString> > mRows;

public:
    PiecewiseFunctionType(const QVector<QVector<QString> >& data);
    QString toString() const override;
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) const;

    QVector<QVector<QString> > data() const;
};


class DataTypeVisitor
{
public:
    virtual void visit(const ConstantType& type) = 0;
    virtual void visit(const FunctionType& type) = 0;
    virtual void visit(const TableType& type) = 0;
    virtual void visit(const PiecewiseFunctionType& type) = 0;
    virtual void visit(const VariableLinkType& type) = 0;
};

#endif // DATATYPE_H
