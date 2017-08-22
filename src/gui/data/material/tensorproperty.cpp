#include "tensorproperty.h"

TensorProperty::TensorProperty() : NamedEntity()
{
    this->mUnit = "";
}

TensorProperty::TensorProperty(const QString &name, const QString &unit) :
    NamedEntity(name)
{
    this->mUnit = unit;
}

TensorProperty::TensorProperty(const TensorProperty &tp) : NamedEntity(tp)
{
    this->mUnit = tp.mUnit;
    this->mIM = tp.mIM;
    this->mDM = tp.mDM;
    this->mSM = tp.mSM;
    this->mAM = tp.mAM;
}

void TensorProperty::setIsotropicModel(const IsotropicModel& im)
{
    this->mIM = im;
}

void TensorProperty::setDiagonalModel(const DiagonalModel& dm)
{
    this->mDM = dm;
}

void TensorProperty::setSymmetricModel(const SymmetricModel& sm)
{
    this->mSM = sm;
}

void TensorProperty::setAnisotropicModel(const AnisotropicModel& am)
{
    this->mAM = am;
}

QString TensorProperty::unit() const
{
    return this->mUnit;
}

IsotropicModel TensorProperty::isotropicModel() const
{
    return this->mIM;
}

DiagonalModel TensorProperty::diagonalModel() const
{
    return this->mDM;
}

SymmetricModel TensorProperty::symmetricModel() const
{
    return this->mSM;
}

AnisotropicModel TensorProperty::anisotropicModel() const
{
    return this->mAM;
}
