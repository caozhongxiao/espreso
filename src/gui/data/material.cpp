#include "material.h"

MaterialPropertyMatrix::~MaterialPropertyMatrix()
{
    delete this->kxx;
}

MaterialPropertyMatrix* BasicMatrix::copy() const
{
    BasicMatrix* matrix = new BasicMatrix();
    matrix->kxx = this->kxx->copy();
    return matrix;
}

MaterialPropertyMatrix* IsotropicMatrix::copy() const
{
    IsotropicMatrix* matrix = new IsotropicMatrix();
    matrix->kxx = this->kxx->copy();
    return matrix;
}

DiagonalMatrix::~DiagonalMatrix()
{
    delete this->kyy;
    delete this->kzz;
}

MaterialPropertyMatrix* DiagonalMatrix::copy() const
{
    DiagonalMatrix* matrix = new DiagonalMatrix();
    matrix->kxx = this->kxx->copy();
    matrix->kyy = this->kyy->copy();
    matrix->kzz = this->kzz->copy();
    return matrix;
}

SymmetricMatrix::~SymmetricMatrix()
{
    delete this->kyy;
    delete this->kzz;
    delete this->kxy;
    delete this->kxz;
    delete this->kyz;
}

MaterialPropertyMatrix* SymmetricMatrix::copy() const
{
    SymmetricMatrix* matrix = new SymmetricMatrix();
    matrix->kxx = this->kxx->copy();
    matrix->kyy = this->kyy->copy();
    matrix->kzz = this->kzz->copy();
    matrix->kxy = this->kxy->copy();
    matrix->kxz = this->kxz->copy();
    matrix->kyz = this->kyz->copy();
    return matrix;
}

AnisotropicMatrix::~AnisotropicMatrix()
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

MaterialPropertyMatrix* AnisotropicMatrix::copy() const
{
    AnisotropicMatrix* matrix = new AnisotropicMatrix();
    matrix->kxx = this->kxx->copy();
    matrix->kyy = this->kyy->copy();
    matrix->kzz = this->kzz->copy();
    matrix->kxy = this->kxy->copy();
    matrix->kxz = this->kxz->copy();
    matrix->kyz = this->kyz->copy();
    matrix->kyx = this->kyx->copy();
    matrix->kzx = this->kzx->copy();
    matrix->kzy = this->kzy->copy();
    return matrix;
}

MaterialProperty::MaterialProperty()
{
    this->mName = "";
    this->mUnit = "";
    this->mValue = nullptr;
}

MaterialProperty::MaterialProperty(const QString& name, const QString& unit, const MaterialPropertyMatrix* value) :
    NamedEntity(name)
{
    this->mUnit = unit;
    this->mValue = value->copy();
}

MaterialProperty::MaterialProperty(const MaterialProperty& mp) :
    NamedEntity(mp)
{
    this->mUnit = mp.mUnit;
    this->mValue = mp.mValue->copy();
}

MaterialProperty::~MaterialProperty()
{
    delete this->mValue;
}

const QString& MaterialProperty::unit() const
{
    return this->mUnit;
}


const MaterialPropertyMatrix* MaterialProperty::value() const
{
    return this->mValue;
}


BasicProperty::BasicProperty() : MaterialProperty() {}

BasicProperty::BasicProperty(const QString& name, const QString& unit, const BasicMatrix* value) :
    MaterialProperty(name, unit, value) {}

BasicProperty::BasicProperty(const BasicProperty& bp) : MaterialProperty(bp) {}

void BasicProperty::accept(MaterialPropertyVisitor* visitor) const
{
    visitor->visit(*this);
}

MaterialProperty* BasicProperty::copy() const
{
    return new BasicProperty(*this);
}


IsotropicProperty::IsotropicProperty() : MaterialProperty() {}

IsotropicProperty::IsotropicProperty(const QString& name, const QString& unit, const IsotropicMatrix* value) :
    MaterialProperty(name, unit, value) {}

IsotropicProperty::IsotropicProperty(const IsotropicProperty& ip) : MaterialProperty(ip) {}

void IsotropicProperty::accept(MaterialPropertyVisitor* visitor) const
{
    visitor->visit(*this);
}

MaterialProperty* IsotropicProperty::copy() const
{
    return new IsotropicProperty(*this);
}


DiagonalProperty::DiagonalProperty() : MaterialProperty() {}

DiagonalProperty::DiagonalProperty(const QString& name, const QString& unit, const DiagonalMatrix* value) :
    MaterialProperty(name, unit, value) {}

DiagonalProperty::DiagonalProperty(const DiagonalProperty& dp) : MaterialProperty(dp) {}

void DiagonalProperty::accept(MaterialPropertyVisitor* visitor) const
{
    visitor->visit(*this);
}

MaterialProperty* DiagonalProperty::copy() const
{
    return new DiagonalProperty(*this);
}

SymmetricProperty::SymmetricProperty() : MaterialProperty() {}

SymmetricProperty::SymmetricProperty(const QString& name, const QString& unit, const SymmetricMatrix* value) :
    MaterialProperty(name, unit, value) {}

SymmetricProperty::SymmetricProperty(const SymmetricProperty& sp) : MaterialProperty(sp) {}

void SymmetricProperty::accept(MaterialPropertyVisitor* visitor) const
{
    visitor->visit(*this);
}

MaterialProperty* SymmetricProperty::copy() const
{
    return new SymmetricProperty(*this);
}


AnisotropicProperty::AnisotropicProperty() : MaterialProperty() {}

AnisotropicProperty::AnisotropicProperty(const QString& name, const QString& unit, const AnisotropicMatrix* value) :
    MaterialProperty(name, unit, value) {}

AnisotropicProperty::AnisotropicProperty(const AnisotropicProperty& ap) : MaterialProperty(ap) {}

void AnisotropicProperty::accept(MaterialPropertyVisitor* visitor) const
{
    visitor->visit(*this);
}

MaterialProperty* AnisotropicProperty::copy() const
{
    return new AnisotropicProperty(*this);
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
    foreach (MaterialProperty* mp, m.mProperties) {
        this->mProperties.append(mp->copy());
    }
}

Material::~Material()
{
    foreach (MaterialProperty* mp, this->mProperties) {
        delete mp;
    }
    this->mProperties.clear();
}

QString Material::description() const
{
    return this->mDescription;
}

void Material::appendProperty(MaterialProperty* const property)
{
    this->mProperties.append(property->copy());
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

const MaterialProperty* Material::property(int index) const
{
    if (index >= this->mProperties.count())
        return nullptr;
    return this->mProperties.at(index);
}

QString Material::toString() const
{
    return this->name();
}
