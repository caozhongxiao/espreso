#include "scalarproperty.h"

ScalarProperty::ScalarProperty() : NamedEntity()
{
    this->mUnit = "";
    this->mAbbreviation = "";
}

ScalarProperty::ScalarProperty(const QString& name, const QString& unit, const QString& abbreviation,
                               DataType* data) :
    NamedEntity(name)
{
    this->mUnit = unit;
    this->mAbbreviation = abbreviation;
    this->mData = data;
}

ScalarProperty::ScalarProperty(const ScalarProperty& sp) : NamedEntity(sp)
{
    this->mUnit = sp.mUnit;
    this->mAbbreviation = sp.mAbbreviation;
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

QString ScalarProperty::abbreviation() const
{
    return this->mAbbreviation;
}
