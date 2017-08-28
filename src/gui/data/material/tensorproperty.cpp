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
    this->mModels = tp.mModels;
    this->mActiveModel = tp.mActiveModel;
}

int TensorProperty::activeModel()
{
    return this->mActiveModel;
}

QString TensorProperty::unit() const
{
    return this->mUnit;
}

void TensorProperty::appendModel(const TensorPropertyModel& model)
{
    this->mModels.append(model);
}

const TensorPropertyModel& TensorProperty::model(int index) const
{
    return this->mModels.at(index);
}

auto TensorProperty::modelBegin()
{
    return this->mModels.begin();
}

auto TensorProperty::modelEnd()
{
    return this->mModels.end();
}
