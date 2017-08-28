#include "tensorproperty.h"

TensorProperty::TensorProperty() : NamedEntity()
{
    this->mActiveModel = 0;
}

TensorProperty::TensorProperty(const QString &name) :
    NamedEntity(name)
{
    this->mActiveModel = 0;
}

TensorProperty::TensorProperty(const TensorProperty &tp) : NamedEntity(tp)
{
    this->mModels = tp.mModels;
    this->mActiveModel = tp.mActiveModel;
}

int TensorProperty::activeModel()
{
    return this->mActiveModel;
}

void TensorProperty::setActiveModel(int model)
{
    this->mActiveModel = model;
}

void TensorProperty::appendModel(const TensorPropertyModel& model)
{
    this->mModels.append(model);
}

const TensorPropertyModel& TensorProperty::model(int index) const
{
    return this->mModels.at(index);
}

QVector<TensorPropertyModel>::Iterator TensorProperty::modelBegin()
{
    return this->mModels.begin();
}

QVector<TensorPropertyModel>::Iterator TensorProperty::modelEnd()
{
    return this->mModels.end();
}
