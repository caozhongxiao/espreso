#include "namedentity.h"

NamedEntity::NamedEntity()
{
    this->mName = "";
}

NamedEntity::NamedEntity(const QString& name)
{
    this->mName = name;
}

QString NamedEntity::name() const
{
    return this->mName;
}
