#include "material.h"

MaterialProperty::MaterialProperty() : NamedEntity()
{
    this->mUnit = "";
}

MaterialProperty::MaterialProperty(const QString& name, const QString& unit) :
    NamedEntity(name)
{
    this->mUnit = unit;
}

MaterialProperty::MaterialProperty(const MaterialProperty& mp) :
    NamedEntity(mp)
{
    this->mUnit = mp.mUnit;
}

const QString& MaterialProperty::unit() const
{
    return this->mUnit;
}


BasicProperty::BasicProperty() : MaterialProperty() {}

BasicProperty::BasicProperty(const QString& name, const QString& unit, const BasicModel& model) :
    MaterialProperty(name, unit)
{
    this->mModel = model;
}

BasicProperty::BasicProperty(const BasicProperty& bp) : MaterialProperty(bp)
{
    this->mModel = bp.mModel;
}

BasicModel& BasicProperty::model()
{
    return this->mModel;
}

void BasicProperty::setModel(const BasicModel& model)
{
    this->mModel = model;
}

QVector<DataType*> BasicProperty::modelData()
{
    QVector<DataType*> model;
    model.append(this->mModel.kxx);
    for (int i = 0; i < 8; ++i)
        model.append(new DummyType());
    return model;
}

int BasicProperty::setModelData(const QVector<DataType*>& data)
{
    if (data.size() != 1)
        return -1;
    mModel.kxx = data.at(0);
    return 1;
}

void BasicProperty::accept(MaterialPropertyVisitor* visitor)
{
    visitor->visit(*this);
}


IsotropicProperty::IsotropicProperty() : MaterialProperty() {}

IsotropicProperty::IsotropicProperty(const QString& name, const QString& unit, const IsotropicModel& model) :
    MaterialProperty(name, unit)
{
    this->mModel = model;
}

IsotropicProperty::IsotropicProperty(const IsotropicProperty& p) : MaterialProperty(p)
{
    this->mModel = p.mModel;
}

IsotropicModel& IsotropicProperty::model()
{
    return this->mModel;
}

void IsotropicProperty::setModel(const IsotropicModel& model)
{
    this->mModel = model;
}

QVector<DataType*> IsotropicProperty::modelData()
{
    QVector<DataType*> model;
    model.append(this->mModel.kxx);
    for (int i = 0; i < 8; ++i)
        model.append(new DummyType());
    return model;
}

int IsotropicProperty::setModelData(const QVector<DataType*>& data)
{
    if (data.size() != 1)
        return -1;
    mModel.kxx = data.at(0);
    return 1;
}

void IsotropicProperty::accept(MaterialPropertyVisitor* visitor)
{
    visitor->visit(*this);
}

DiagonalProperty::DiagonalProperty() : MaterialProperty() {}

DiagonalProperty::DiagonalProperty(const QString& name, const QString& unit, const DiagonalModel& model) :
    MaterialProperty(name, unit)
{
    this->mModel = model;
}

DiagonalProperty::DiagonalProperty(const DiagonalProperty& p) : MaterialProperty(p)
{
    this->mModel = p.mModel;
}

DiagonalModel& DiagonalProperty::model()
{
    return this->mModel;
}

void DiagonalProperty::setModel(const DiagonalModel& model)
{
    this->mModel = model;
}

QVector<DataType*> DiagonalProperty::modelData()
{
    QVector<DataType*> model;
    model.append(this->mModel.kxx);
    for (int i = 0; i < 3; ++i)
        model.append(new DummyType());
    model.append(this->mModel.kyy);
    for (int i = 0; i < 3; ++i)
        model.append(new DummyType());
    model.append(this->mModel.kzz);
    return model;
}

int DiagonalProperty::setModelData(const QVector<DataType*>& data)
{
    if (data.size() != 3)
        return -1;
    mModel.kxx = data.at(0);
    mModel.kyy = data.at(1);
    mModel.kzz = data.at(2);
    return 3;
}

void DiagonalProperty::accept(MaterialPropertyVisitor* visitor)
{
    visitor->visit(*this);
}


SymmetricProperty::SymmetricProperty() : MaterialProperty() {}

SymmetricProperty::SymmetricProperty(const QString& name, const QString& unit, const SymmetricModel& model) :
    MaterialProperty(name, unit)
{
    this->mModel = model;
}

SymmetricProperty::SymmetricProperty(const SymmetricProperty& p) : MaterialProperty(p)
{
    this->mModel = p.mModel;
}

SymmetricModel& SymmetricProperty::model()
{
    return this->mModel;
}

void SymmetricProperty::setModel(const SymmetricModel& model)
{
    this->mModel = model;
}

QVector<DataType*> SymmetricProperty::modelData()
{
    QVector<DataType*> model;
    model.append(this->mModel.kxx);
    model.append(this->mModel.kxy);
    model.append(this->mModel.kxz);
    model.append(new DummyType());
    model.append(this->mModel.kyy);
    model.append(this->mModel.kyz);
    for (int i = 0; i < 2; ++i)
        model.append(new DummyType());
    model.append(this->mModel.kzz);
    return model;
}

int SymmetricProperty::setModelData(const QVector<DataType*>& data)
{
    if (data.size() != 6)
        return -1;
    mModel.kxx = data.at(0);
    mModel.kyy = data.at(1);
    mModel.kzz = data.at(2);
    mModel.kxy = data.at(3);
    mModel.kxz = data.at(4);
    mModel.kyz = data.at(5);
    return 6;
}

void SymmetricProperty::accept(MaterialPropertyVisitor* visitor)
{
    visitor->visit(*this);
}


AnisotropicProperty::AnisotropicProperty() : MaterialProperty() {}

AnisotropicProperty::AnisotropicProperty(const QString& name, const QString& unit, const AnisotropicModel& model) :
    MaterialProperty(name, unit)
{
    this->mModel = model;
}

AnisotropicProperty::AnisotropicProperty(const AnisotropicProperty& p) : MaterialProperty(p)
{
    this->mModel = p.mModel;
}

AnisotropicModel& AnisotropicProperty::model()
{
    return this->mModel;
}

void AnisotropicProperty::setModel(const AnisotropicModel& model)
{
    this->mModel = model;
}

QVector<DataType*> AnisotropicProperty::modelData()
{
    QVector<DataType*> model;
    model.append(this->mModel.kxx);
    model.append(this->mModel.kxy);
    model.append(this->mModel.kxz);
    model.append(this->mModel.kyx);
    model.append(this->mModel.kyy);
    model.append(this->mModel.kyz);
    model.append(this->mModel.kzx);
    model.append(this->mModel.kzy);
    model.append(this->mModel.kzz);
    return model;
}

int AnisotropicProperty::setModelData(const QVector<DataType*>& data)
{
    if (data.size() != 9)
        return -1;
    mModel.kxx = data.at(0);
    mModel.kyy = data.at(1);
    mModel.kzz = data.at(2);
    mModel.kxy = data.at(3);
    mModel.kxz = data.at(4);
    mModel.kyz = data.at(5);
    mModel.kyx = data.at(6);
    mModel.kzx = data.at(7);
    mModel.kzy = data.at(8);
    return 9;
}

void AnisotropicProperty::accept(MaterialPropertyVisitor* visitor)
{
    visitor->visit(*this);
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
    this->mProperties = m.mProperties;
}

QString Material::description() const
{
    return this->mDescription;
}

void Material::appendProperty(MaterialProperty* property)
{
    this->mProperties.append(property);
}


void Material::removeProperty(int index)
{
    MaterialProperty* tmp = this->mProperties.at(index);
    this->mProperties.remove(index);
    delete tmp;
}

int Material::modifyProperty(int index, MaterialProperty* property)
{
    if (this->mProperties.count() <= index)
        return -1;
    MaterialProperty* tmp = this->mProperties.at(index);
    this->mProperties.replace(index, property);
    delete tmp;
    return 1;
}

MaterialProperty* Material::property(int index) const
{
    if (index >= this->mProperties.count())
        return nullptr;
    return this->mProperties.at(index);
}

QString Material::toString() const
{
    return this->name();
}
