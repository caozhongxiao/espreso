#ifndef MATERIALPROPERTY_H
#define MATERIALPROPERTY_H

#include "../namedentity.h"
#include "../datatype.h"

class MaterialProperty : public NamedEntity
{

public:
    MaterialProperty();
    MaterialProperty(const QString& name, const QString& unit, const QString& symbol,
                   DataType* data);
    MaterialProperty(const MaterialProperty& sp);
    virtual ~MaterialProperty();

    void setData(DataType* data);
    DataType* data() const;

    QString unit() const;
    QString symbol() const;

protected:
    QString mUnit;
    QString mSymbol;
    DataType* mData;
};

#endif // MATERIALPROPERTY_H
