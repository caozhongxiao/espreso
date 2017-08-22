#include "material.h"

Material::Material()
{
    this->mName = "";
    this->mDescription = "";
}

Material::Material(const QString& name, const QString& description) : NamedEntity(name)
{
    this->mDescription = description;
}

Material::Material(const Material& m) : NamedEntity(m)
{
    this->mDescription = m.mDescription;
    this->mScalars = m.mScalars;
    this->mTensors = m.mTensors;
}

QString Material::description() const
{
    return this->mDescription;
}

int Material::append(const ScalarProperty& property)
{
    this->mScalars.append(property);
    return this->mScalars.size();
}

int Material::append(const TensorProperty& property)
{
    this->mTensors.append(property);
    return this->mTensors.size();
}

void Material::removeScalar(int index)
{
    this->mScalars.remove(index);
}

void Material::removeTensor(int index)
{
    this->mTensors.remove(index);
}

ScalarProperty& Material::scalar(int index) const
{
    this->mScalars.at(index);
}

TensorProperty& Material::tensor(int index) const
{
    this->mTensors.at(index);
}
