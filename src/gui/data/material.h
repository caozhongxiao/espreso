#ifndef MATERIAL_H
#define MATERIAL_H

#include <QString>
#include <QVector>
#include "datatype.h"
#include "namedentity.h"

class MaterialPropertyModelVisitor;

struct MaterialPropertyModel
{
    DataType* kxx;
    virtual ~MaterialPropertyModel();
    virtual MaterialPropertyModel* copy() const = 0;
    virtual void accept(MaterialPropertyModelVisitor* visitor) const = 0;
};

struct BasicModel final : MaterialPropertyModel
{
    virtual MaterialPropertyModel* copy() const override;
    virtual void accept(MaterialPropertyModelVisitor* visitor) const override;
};

struct IsotropicModel final : MaterialPropertyModel
{
    virtual MaterialPropertyModel* copy() const override;
    virtual void accept(MaterialPropertyModelVisitor* visitor) const override;
};

struct DiagonalModel final : MaterialPropertyModel
{
    DataType* kyy;
    DataType* kzz;
    ~DiagonalModel();
    virtual MaterialPropertyModel* copy() const override;
    virtual void accept(MaterialPropertyModelVisitor* visitor) const override;
};

struct SymmetricModel final : MaterialPropertyModel
{
    DataType* kyy;
    DataType* kzz;
    DataType* kxy;
    DataType* kxz;
    DataType* kyz;
    ~SymmetricModel();
    virtual MaterialPropertyModel* copy() const override;
    virtual void accept(MaterialPropertyModelVisitor* visitor) const override;
};

struct AnisotropicModel final : MaterialPropertyModel
{
    DataType* kyy;
    DataType* kzz;
    DataType* kxy;
    DataType* kxz;
    DataType* kyz;
    DataType* kyx;
    DataType* kzx;
    DataType* kzy;
    ~AnisotropicModel();
    virtual MaterialPropertyModel* copy() const override;
    virtual void accept(MaterialPropertyModelVisitor* visitor) const override;
};


class MaterialPropertyModelVisitor
{
public:
    virtual void visit(const BasicModel& property) = 0;
    virtual void visit(const IsotropicModel& property) = 0;
    virtual void visit(const DiagonalModel& property) = 0;
    virtual void visit(const SymmetricModel& property) = 0;
    virtual void visit(const AnisotropicModel& property) = 0;
};

class MaterialProperty : public NamedEntity
{
protected:
    QString mUnit;
    MaterialPropertyModel* mModel;

public:
    MaterialProperty();
    MaterialProperty(const QString& name, const QString& unit, MaterialPropertyModel* const model);
    MaterialProperty(const MaterialProperty&);
    virtual ~MaterialProperty();

    const QString& unit() const;
    void replaceModel(MaterialPropertyModel* model);
    MaterialPropertyModel* model() const;
};


class Material : public NamedEntity
{
private:
    QString mDescription;
    QVector<MaterialProperty> mProperties;

public:
    Material();
    Material(const QString& name, const QString& description);
    Material(const Material&);
    ~Material();

    QString description() const;
    QString toString() const;
    void appendProperty(const MaterialProperty& property);
    void removeProperty(int index);
    int modifyProperty(int index, const MaterialProperty& property);
    const MaterialProperty& property(int index) const;
};

#endif // MATERIAL_H
