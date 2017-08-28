#include "materialproperty.h"

MaterialProperty::MaterialProperty() : NamedEntity()
{
    this->mUnit = "";
    this->mSymbol = "";
    this->mData = new DummyType();
}

MaterialProperty::MaterialProperty(const QString& name, const QString& unit, const QString& symbol,
                               DataType* data) :
    NamedEntity(name)
{
    this->mUnit = unit;
    this->mSymbol = symbol;
    this->mData = data;
}

MaterialProperty::MaterialProperty(const MaterialProperty& sp) : NamedEntity(sp)
{
    this->mUnit = sp.mUnit;
    this->mSymbol = sp.mSymbol;
    this->mData = sp.mData->copy();
}

MaterialProperty::~MaterialProperty()
{
    delete this->mData;
}

void MaterialProperty::setData(DataType* data)
{
    if (data == this->mData)
    {
        this->mData = data;
    }
    else
    {
        DataType* tmp = this->mData;
        this->mData = data;
        delete tmp;
    }
}

DataType* MaterialProperty::data() const
{
    return this->mData;
}

QString MaterialProperty::unit() const
{
    return this->mUnit;
}

QString MaterialProperty::symbol() const
{
    return this->mSymbol;
}
