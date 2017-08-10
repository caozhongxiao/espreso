#ifndef NAMEDENTITY_H
#define NAMEDENTITY_H

#include <QString>

class NamedEntity
{
protected:
    QString mName;

    NamedEntity();
    NamedEntity(const QString&);
    NamedEntity(const NamedEntity&);

public:
    const QString& name() const;
};

#endif // NAMEDENTITY_H
