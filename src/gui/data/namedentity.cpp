#include "namedentity.h"

NamedEntity::NamedEntity()
{
    this->mName = "";
}

NamedEntity::NamedEntity(const QString& name)
{
    this->mName = name;
}

NamedEntity::NamedEntity(const NamedEntity& ne)
{
    this->mName = ne.mName;
}

const QString& NamedEntity::name() const
{
    return this->mName;
}
