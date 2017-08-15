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
    new DummyType();
}

void DummyType::accept(DataTypeVisitor *visitor)
{

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


DataType* VariableLinkType::copy() const
{
    return new VariableLinkType(this->data);
}

void VariableLinkType::accept(DataTypeVisitor* visitor)
{
    visitor->visit(*this);
}


TableType::TableType(const QVector<QPair<QString, QString> >& data)
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

QVector<QPair<QString, QString> > TableType::data() const
{
    return this->mRows;
}


PiecewiseFunctionType::PiecewiseFunctionType(const QVector<QVector<QString> >& data)
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

QVector<QVector<QString> > PiecewiseFunctionType::data() const
{
    return this->mRows;
}
