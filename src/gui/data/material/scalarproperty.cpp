#include "scalarproperty.h"

ScalarProperty::ScalarProperty() : NamedEntity()
{
    this->mUnit = "";
    this->mSymbol = "";
}

ScalarProperty::ScalarProperty(const QString& name, const QString& unit, const QString& symbol,
                               DataType* data) :
    NamedEntity(name)
{
    this->mUnit = unit;
    this->mSymbol = symbol;
    this->mData = data;
}

ScalarProperty::ScalarProperty(const ScalarProperty& sp) : NamedEntity(sp)
{
    this->mUnit = sp.mUnit;
    this->mSymbol = sp.mSymbol;
    this->mData = sp.mData;
}

void ScalarProperty::setData(DataType* data)
{
    this->mData = data;
}

DataType* ScalarProperty::data() const
{
    return this->mData;
}

QString ScalarProperty::unit() const
{
    return this->mUnit;
}

QString ScalarProperty::symbol() const
{
    return this->mSymbol;
}
