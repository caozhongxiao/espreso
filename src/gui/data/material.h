#ifndef MATERIAL_H
#define MATERIAL_H

#include <QString>
#include <QVector>
#include "datatype.h"
#include "namedentity.h"


struct MaterialPropertyMatrix
{
    DataType* kxx;
};

struct BasicMatrix final : MaterialPropertyMatrix
{
};

struct IsotropicMatrix final : MaterialPropertyMatrix
{
};

struct DiagonalMatrix final : MaterialPropertyMatrix
{
    DataType* kyy;
    DataType* kzz;
};

struct SymmetricMatrix final : MaterialPropertyMatrix
{
    DataType* kyy;
    DataType* kzz;
    DataType* kxy;
    DataType* kxz;
    DataType* kyz;
};

struct AnisotropicMatrix final : MaterialPropertyMatrix
{
    DataType* kyy;
    DataType* kzz;
    DataType* kxy;
    DataType* kxz;
    DataType* kyz;
    DataType* kyx;
    DataType* kzx;
    DataType* kzy;
};


class MaterialProperty : public NamedEntity
{
protected:
    QString mName;
    QString mUnit;
    MaterialPropertyMatrix* mValue;

    MaterialProperty();
    MaterialProperty(const QString& name, const QString& unit, MaterialPropertyMatrix* const value);

public:
    const QString& name() const;
    const QString& unit() const;
    MaterialPropertyMatrix* value() const;
};

class BasicProperty final : public MaterialProperty
{
public:
    BasicProperty();
    BasicProperty(const QString& name, const QString& unit, BasicMatrix* const value);
};

class IsotropicProperty final : public MaterialProperty
{
public:
    IsotropicProperty();
    IsotropicProperty(const QString& name, const QString& unit, IsotropicMatrix* const value);
};

class DiagonalProperty final : public MaterialProperty
{
public:
    DiagonalProperty();
    DiagonalProperty(const QString& name, const QString& unit, DiagonalMatrix* const value);
};

class SymmetricProperty final : public MaterialProperty
{
public:
    SymmetricProperty();
    SymmetricProperty(const QString& name, const QString& unit, SymmetricMatrix* const value);
};

class AnisotropicProperty final : public MaterialProperty
{
public:
    AnisotropicProperty();
    AnisotropicProperty(const QString& name, const QString& unit, AnisotropicMatrix* const value);
};


class Material : public NamedEntity
{
private:
    QString mName;
    QString mDescription;
    QVector<MaterialProperty*> mProperties;

public:
    Material();
    Material(const QString& name, const QString& description);

    QString description() const;
    QString name() const;
    void appendProperty(MaterialProperty* const property);
    void removeProperty(int index);
    int modifyProperty(int index, MaterialProperty* const property);
    MaterialProperty* property(int index) const;
};

#endif // MATERIAL_H
