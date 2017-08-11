#include "material.h"

MaterialPropertyModel::~MaterialPropertyModel()
{
    delete this->kxx;
}

MaterialPropertyModel* BasicModel::copy() const
{
    BasicModel* Model = new BasicModel();
    Model->kxx = this->kxx->copy();
    return Model;
}

void BasicModel::accept(MaterialPropertyModelVisitor *visitor) const
{
    visitor->visit(*this);
}

MaterialPropertyModel* IsotropicModel::copy() const
{
    IsotropicModel* Model = new IsotropicModel();
    Model->kxx = this->kxx->copy();
    return Model;
}

void IsotropicModel::accept(MaterialPropertyModelVisitor *visitor) const
{
    visitor->visit(*this);
}

DiagonalModel::~DiagonalModel()
{
    delete this->kyy;
    delete this->kzz;
}

MaterialPropertyModel* DiagonalModel::copy() const
{
    DiagonalModel* Model = new DiagonalModel();
    Model->kxx = this->kxx->copy();
    Model->kyy = this->kyy->copy();
    Model->kzz = this->kzz->copy();
    return Model;
}

void DiagonalModel::accept(MaterialPropertyModelVisitor *visitor) const
{
    visitor->visit(*this);
}

SymmetricModel::~SymmetricModel()
{
    delete this->kyy;
    delete this->kzz;
    delete this->kxy;
    delete this->kxz;
    delete this->kyz;
}

MaterialPropertyModel* SymmetricModel::copy() const
{
    SymmetricModel* Model = new SymmetricModel();
    Model->kxx = this->kxx->copy();
    Model->kyy = this->kyy->copy();
    Model->kzz = this->kzz->copy();
    Model->kxy = this->kxy->copy();
    Model->kxz = this->kxz->copy();
    Model->kyz = this->kyz->copy();
    return Model;
}

void SymmetricModel::accept(MaterialPropertyModelVisitor *visitor) const
{
    visitor->visit(*this);
}

AnisotropicModel::~AnisotropicModel()
{
    delete this->kyy;
    delete this->kzz;
    delete this->kxy;
    delete this->kxz;
    delete this->kyz;
    delete this->kyx;
    delete this->kzx;
    delete this->kzy;
}

MaterialPropertyModel* AnisotropicModel::copy() const
{
    AnisotropicModel* Model = new AnisotropicModel();
    Model->kxx = this->kxx->copy();
    Model->kyy = this->kyy->copy();
    Model->kzz = this->kzz->copy();
    Model->kxy = this->kxy->copy();
    Model->kxz = this->kxz->copy();
    Model->kyz = this->kyz->copy();
    Model->kyx = this->kyx->copy();
    Model->kzx = this->kzx->copy();
    Model->kzy = this->kzy->copy();
    return Model;
}

void AnisotropicModel::accept(MaterialPropertyModelVisitor *visitor) const
{
    visitor->visit(*this);
}


MaterialProperty::MaterialProperty() : NamedEntity()
{
    this->mUnit = "";
    this->mModel = nullptr;
}

MaterialProperty::MaterialProperty(const QString& name, const QString& unit, MaterialPropertyModel* const model) :
    NamedEntity(name)
{
    this->mUnit = unit;
    this->mModel = model;
}

MaterialProperty::MaterialProperty(const MaterialProperty& mp) :
    NamedEntity(mp)
{
    this->mUnit = mp.mUnit;
    this->mModel = mp.mModel->copy();
}

MaterialProperty::~MaterialProperty()
{
    delete this->mModel;
}

const QString& MaterialProperty::unit() const
{
    return this->mUnit;
}


MaterialPropertyModel* MaterialProperty::model() const
{
    return this->mModel;
}


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
    foreach (MaterialProperty mp, m.mProperties) {
        this->mProperties.append(MaterialProperty(mp));
    }
}

Material::~Material()
{
    this->mProperties.clear();
}

QString Material::description() const
{
    return this->mDescription;
}

void Material::appendProperty(const MaterialProperty& property)
{
    this->mProperties.append(property);
}


void Material::removeProperty(int index)
{
    this->mProperties.remove(index);
}

int Material::modifyProperty(int index, const MaterialProperty& property)
{
    if (this->mProperties.count() <= index)
        return -1;
    this->mProperties.replace(index, property);
    return 1;
}

const MaterialProperty& Material::property(int index) const
{
    if (index >= this->mProperties.count())
        return MaterialProperty();
    return this->mProperties.at(index);
}

QString Material::toString() const
{
    return this->name();
}
