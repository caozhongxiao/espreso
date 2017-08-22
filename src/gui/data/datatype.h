#ifndef DATATYPE_H
#define DATATYPE_H

#include <QString>
#include <QVariant>
#include <QVector>
#include <QList>
#include <QObject>
#include <QStringList>

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
    virtual void accept(DataTypeVisitor* visitor) = 0;
};

class StringType : public DataType
{
protected:
    QString data;
    StringType(const QString& data);

public:
    QString toString() const override;
    virtual DataType* copy() const = 0;
    virtual void accept(DataTypeVisitor* visitor) = 0;
};

class DummyType : public DataType
{
public:
    DummyType() {}
    QString toString() const override;
    DataType* copy() const override;
    void accept(DataTypeVisitor *visitor) override;
};

class ConstantType : public StringType
{
public:
    ConstantType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) override;
};

class FunctionType : public StringType
{
public:
    FunctionType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) override;
};

class ExpressionType : public StringType
{
public:
    ExpressionType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) override;
};

class VariableLinkType : public StringType
{
public:
    VariableLinkType(const QString& data) : StringType(data) {}
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) override;
};

class TableType : public DataType
{
protected:
    QList<QList<QString> > mRows;

public:
    TableType(const QList<QList<QString> >& data);
    QString toString() const override;
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) override;

    const QList<QList<QString> >& data() const;
    static QStringList headlines();
};

class PiecewiseFunctionType : public DataType
{
protected:
    QList<QList<QString> > mRows;

public:
    PiecewiseFunctionType(const QList<QList<QString> >& data);
    QString toString() const override;
    DataType* copy() const override;
    void accept(DataTypeVisitor* visitor) override;

    const QList<QList<QString> >& data() const;
    static QStringList headlines();
};


class DataTypeVisitor
{
public:
    virtual void visit(ConstantType& type) = 0;
    virtual void visit(FunctionType& type) = 0;
    virtual void visit(ExpressionType& type) = 0;
    virtual void visit(TableType& type) = 0;
    virtual void visit(PiecewiseFunctionType& type) = 0;
    virtual void visit(VariableLinkType& type) = 0;

    static QStringList types();
};

#endif // DATATYPE_H
