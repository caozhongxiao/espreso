#include "material.h"

MaterialProperty::MaterialProperty()
{
    this->mName = "";
    this->mUnit = "";
    this->mValue = nullptr;
}

MaterialProperty::MaterialProperty(const QString& name, const QString& unit, MaterialPropertyMatrix* const value)
{
    this->mName = name;
    this->mUnit = unit;
    this->mValue = value;
}

const QString& MaterialProperty::name() const
{
    return this->mName;
}

const QString& MaterialProperty::unit() const
{
    return this->mUnit;
}


MaterialPropertyMatrix* MaterialProperty::value() const
{
    return this->mValue;
}


BasicProperty::BasicProperty() : MaterialProperty() {}

BasicProperty::BasicProperty(const QString& name, const QString& unit, BasicMatrix* const value) :
    MaterialProperty(name, unit, value) {}


IsotropicProperty::IsotropicProperty() : MaterialProperty() {}

IsotropicProperty::IsotropicProperty(const QString& name, const QString& unit, IsotropicMatrix* const value) :
    MaterialProperty(name, unit, value) {}


DiagonalProperty::DiagonalProperty() : MaterialProperty() {}

DiagonalProperty::DiagonalProperty(const QString& name, const QString& unit, DiagonalMatrix* const value) :
    MaterialProperty(name, unit, value) {}


SymmetricProperty::SymmetricProperty() : MaterialProperty() {}

SymmetricProperty::SymmetricProperty(const QString& name, const QString& unit, SymmetricMatrix* const value) :
    MaterialProperty(name, unit, value) {}


AnisotropicProperty::AnisotropicProperty() : MaterialProperty() {}

AnisotropicProperty::AnisotropicProperty(const QString& name, const QString& unit, AnisotropicMatrix* const value) :
    MaterialProperty(name, unit, value) {}


Material::Material()
{
    this->mName = "";
    this->mDescription = "";
}

Material::Material(const QString& name, const QString& description)
{
    this->mName = name;
    this->mDescription = description;
}

QString Material::description() const
{
    return this->mDescription;
}

QString Material::name() const
{
    return this->mName;
}

void Material::appendProperty(MaterialProperty* const property)
{
    this->mProperties.append(property);
}


void Material::removeProperty(int index)
{
    this->mProperties.remove(index);
}

int Material::modifyProperty(int index, MaterialProperty* const property)
{
    if (this->mProperties.count() <= index)
        return -1;
    this->mProperties.replace(index, property);
    return 1;
}

MaterialProperty* Material::property(int index) const
{
    if (index >= this->mProperties.count())
        return nullptr;
    return this->mProperties.at(index);
}
