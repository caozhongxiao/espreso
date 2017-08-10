#ifndef MATERIAL_H
#define MATERIAL_H

#include <QString>
#include <QVector>
#include "datatype.h"
#include "namedentity.h"


struct MaterialPropertyMatrix
{
    DataType* kxx;
    virtual ~MaterialPropertyMatrix();
    virtual MaterialPropertyMatrix* copy() const = 0;
};

struct BasicMatrix final : MaterialPropertyMatrix
{
    virtual MaterialPropertyMatrix* copy() const override;
};

struct IsotropicMatrix final : MaterialPropertyMatrix
{
    virtual MaterialPropertyMatrix* copy() const override;
};

struct DiagonalMatrix final : MaterialPropertyMatrix
{
    DataType* kyy;
    DataType* kzz;
    ~DiagonalMatrix();
    virtual MaterialPropertyMatrix* copy() const override;
};

struct SymmetricMatrix final : MaterialPropertyMatrix
{
    DataType* kyy;
    DataType* kzz;
    DataType* kxy;
    DataType* kxz;
    DataType* kyz;
    ~SymmetricMatrix();
    virtual MaterialPropertyMatrix* copy() const override;
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
    ~AnisotropicMatrix();
    virtual MaterialPropertyMatrix* copy() const override;
};


class MaterialPropertyVisitor;


class MaterialProperty : public NamedEntity
{
protected:
    QString mUnit;
    MaterialPropertyMatrix* mValue;

    MaterialProperty();
    MaterialProperty(const QString& name, const QString& unit, const MaterialPropertyMatrix* value);
    MaterialProperty(const MaterialProperty&);

public:
    virtual ~MaterialProperty();
    const QString& unit() const;
    const MaterialPropertyMatrix* value() const;
    virtual void accept(MaterialPropertyVisitor* visitor) const = 0;
    virtual MaterialProperty* copy() const = 0;
};

class BasicProperty final : public MaterialProperty
{
public:
    BasicProperty();
    BasicProperty(const QString& name, const QString& unit, const BasicMatrix* value);
    BasicProperty(const BasicProperty&);
    void accept(MaterialPropertyVisitor* visitor) const override;
    virtual MaterialProperty* copy() const override;
};

class IsotropicProperty final : public MaterialProperty
{
public:
    IsotropicProperty();
    IsotropicProperty(const QString& name, const QString& unit, const IsotropicMatrix* value);
    IsotropicProperty(const IsotropicProperty&);
    void accept(MaterialPropertyVisitor* visitor) const override;
    virtual MaterialProperty* copy() const override;
};

class DiagonalProperty final : public MaterialProperty
{
public:
    DiagonalProperty();
    DiagonalProperty(const QString& name, const QString& unit, const DiagonalMatrix* value);
    DiagonalProperty(const DiagonalProperty&);
    void accept(MaterialPropertyVisitor* visitor) const override;
    virtual MaterialProperty* copy() const override;
};

class SymmetricProperty final : public MaterialProperty
{
public:
    SymmetricProperty();
    SymmetricProperty(const QString& name, const QString& unit, const SymmetricMatrix* value);
    SymmetricProperty(const SymmetricProperty&);
    void accept(MaterialPropertyVisitor* visitor) const override;
    virtual MaterialProperty* copy() const override;
};

class AnisotropicProperty final : public MaterialProperty
{
public:
    AnisotropicProperty();
    AnisotropicProperty(const QString& name, const QString& unit, const AnisotropicMatrix* value);
    AnisotropicProperty(const AnisotropicProperty&);
    void accept(MaterialPropertyVisitor* visitor) const override;
    virtual MaterialProperty* copy() const override;
};

class MaterialPropertyVisitor
{
public:
    virtual void visit(const BasicProperty& property) = 0;
    virtual void visit(const IsotropicProperty& property) = 0;
    virtual void visit(const DiagonalProperty& property) = 0;
    virtual void visit(const SymmetricProperty& property) = 0;
    virtual void visit(const AnisotropicProperty& property) = 0;
};


class Material : public NamedEntity
{
private:
    QString mDescription;
    QVector<MaterialProperty*> mProperties;

public:
    Material();
    Material(const QString& name, const QString& description);
    Material(const Material&);
    ~Material();

    QString description() const;
    QString toString() const;
    void appendProperty(MaterialProperty* const property);
    void removeProperty(int index);
    int modifyProperty(int index, MaterialProperty* const property);
    const MaterialProperty* property(int index) const;
};

#endif // MATERIAL_H
