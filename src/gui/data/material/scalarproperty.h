#ifndef SCALARPROPERTIES_H
#define SCALARPROPERTIES_H

#include "../namedentity.h"
#include "../datatype.h"

class ScalarProperty : public NamedEntity
{
private:
    QString mUnit;
    QString mAbbreviation;
    DataType* mData;

public:
    ScalarProperty();
    ScalarProperty(const QString& name, const QString& unit, const QString& abbreviation,
                   DataType* data);
    ScalarProperty(const ScalarProperty& sp);

    void setData(DataType* data);
    DataType* data() const;

    QString unit() const;
    QString abbreviation() const;
};

#endif // SCALARPROPERTIES_H
