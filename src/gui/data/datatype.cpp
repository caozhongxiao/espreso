#include "datatype.h"

StringType::StringType(const QString& data):DataType()
{
    this->data = data;
}

QString StringType::toString() const
{
    return this->data;
}

QString DummyType::toString() const
{
    return "";
}

DataType* DummyType::copy() const
{
    return new DummyType();
}

void DummyType::accept(DataTypeVisitor *visitor)
{
    visitor->visit(*this);
}


DataType* ConstantType::copy() const
{
    return new ConstantType(this->data);
}

void ConstantType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}


DataType* FunctionType::copy() const
{
    return new FunctionType(this->data);
}

void FunctionType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}


DataType* ExpressionType::copy() const
{
    return new ExpressionType(this->data);
}

void ExpressionType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}


DataType* VariableLinkType::copy() const
{
    return new VariableLinkType(this->data);
}

void VariableLinkType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}


TableType::TableType(const QList<QList<QString> >& data)
{
    this->mRows = data;
}

QString TableType::toString() const
{
    return QObject::tr("Table");
}

DataType* TableType::copy() const
{
    return new TableType(this->mRows);
}

void TableType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}

const QList<QList<QString> >& TableType::data() const
{
    return this->mRows;
}

QStringList TableType::headlines()
{
    QStringList result;
    result << QObject::tr("x") << QObject::tr("f(x)");

    return result;
}


PiecewiseFunctionType::PiecewiseFunctionType(const QList<QList<QString> >& data)
{
    this->mRows = data;
}

QString PiecewiseFunctionType::toString() const
{
    return QObject::tr("Piecewise function");
}

DataType* PiecewiseFunctionType::copy() const
{
    return new PiecewiseFunctionType(this->mRows);
}

void PiecewiseFunctionType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}

const QList<QList<QString> >& PiecewiseFunctionType::data() const
{
    return this->mRows;
}

QStringList PiecewiseFunctionType::headlines()
{
    QStringList result;
    result << QObject::tr("Lower bound")
           << QObject::tr("Upper bound")
           << QObject::tr("f(x)");

    return result;
}


QStringList DataTypeVisitor::types()
{
    QStringList result;
    result << QObject::tr("Expression")
           << QObject::tr("Table")
           << QObject::tr("Piecewise function");

    return result;
}
