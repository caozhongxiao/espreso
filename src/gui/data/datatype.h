#ifndef DATATYPE_H
#define DATATYPE_H

#include <QString>
#include <QVariant>
#include <QVector>
#include <QPair>
#include <QObject>

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

class TableType : public DataType
{
protected:
    QVector<QPair<QString, QString> > mRows;

public:
    TableType(const QVector<QPair<QString, QString> >& data);
    QString toString() const override;
    DataType* copy() const override;
    virtual int type() const override;

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
    virtual int type() const override;

    QVector<QVector<QString> > data() const;
};

#endif // DATATYPE_H
