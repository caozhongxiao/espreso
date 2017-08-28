#ifndef SCALARPROPERTIES_H
#define SCALARPROPERTIES_H

#include "materialproperty.h"

class ScalarProperty : public MaterialProperty
{
public:
    ScalarProperty() : MaterialProperty() {}
    ScalarProperty(const QString& name, const QString& unit, const QString& symbol,
                   DataType* data) :
        MaterialProperty(name, unit, symbol, data) {}
    ScalarProperty(const ScalarProperty& sp) : MaterialProperty(sp) {}
};

#endif // SCALARPROPERTIES_H
