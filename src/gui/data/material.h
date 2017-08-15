#ifndef MATERIAL_H
#define MATERIAL_H

#include <QString>
#include <QVector>
#include "datatype.h"
#include "namedentity.h"

//MATRIX CELL ODRDER - ALPHABETICAL
//DataType* kxx;
//DataType* kxy;
//DataType* kxz;
//DataType* kyx;
//DataType* kyy;
//DataType* kyz;
//DataType* kzx;
//DataType* kzy;
//DataType* kzz;

class MaterialPropertyVisitor;

struct BasicModel
{
    DataType* kxx;
};

struct IsotropicModel
{
    DataType* kxx;
};

struct DiagonalModel
{
    DataType* kxx;
    DataType* kyy;
    DataType* kzz;
};

struct SymmetricModel
{
    DataType* kxx;
    DataType* kyy;
    DataType* kzz;
    DataType* kxy;
    DataType* kxz;
    DataType* kyz;
};

struct AnisotropicModel
{
    DataType* kxx;
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
    QString mUnit;
    MaterialProperty();
    MaterialProperty(const QString& name, const QString& unit);
    MaterialProperty(const MaterialProperty&);

public:
    virtual ~MaterialProperty() {}

    const QString& unit() const;
    virtual QVector<DataType*> modelData() = 0;
    virtual int setModelData(const QVector<DataType*>& data) = 0;

    virtual void accept(MaterialPropertyVisitor* visitor) = 0;
};

class BasicProperty : public MaterialProperty
{
private:
    BasicModel mModel;
public:
    BasicProperty();
    BasicProperty(const QString& name, const QString& unit, const BasicModel& model);
    BasicProperty(const BasicProperty&);

    BasicModel& model();
    void setModel(const BasicModel& model);

    QVector<DataType*> modelData() override;
    int setModelData(const QVector<DataType*>& data) override;

    void accept(MaterialPropertyVisitor* visitor) override;
};

class IsotropicProperty : public MaterialProperty
{
private:
    IsotropicModel mModel;
public:
    IsotropicProperty();
    IsotropicProperty(const QString& name, const QString& unit, const IsotropicModel& model);
    IsotropicProperty(const IsotropicProperty&);

    IsotropicModel& model();
    void setModel(const IsotropicModel& model);

    QVector<DataType*> modelData() override;
    int setModelData(const QVector<DataType*>& data) override;

    void accept(MaterialPropertyVisitor* visitor) override;
};

class DiagonalProperty : public MaterialProperty
{
private:
    DiagonalModel mModel;
public:
    DiagonalProperty();
    DiagonalProperty(const QString& name, const QString& unit, const DiagonalModel& model);
    DiagonalProperty(const DiagonalProperty&);

    DiagonalModel& model();
    void setModel(const DiagonalModel& model);

    QVector<DataType*> modelData() override;
    int setModelData(const QVector<DataType*>& data) override;

    void accept(MaterialPropertyVisitor* visitor) override;
};

class SymmetricProperty : public MaterialProperty
{
private:
    SymmetricModel mModel;
public:
    SymmetricProperty();
    SymmetricProperty(const QString& name, const QString& unit, const SymmetricModel& model);
    SymmetricProperty(const SymmetricProperty&);

    SymmetricModel& model();
    void setModel(const SymmetricModel& model);

    QVector<DataType*> modelData() override;
    int setModelData(const QVector<DataType*>& data) override;

    void accept(MaterialPropertyVisitor* visitor) override;
};

class AnisotropicProperty : public MaterialProperty
{
private:
    AnisotropicModel mModel;
public:
    AnisotropicProperty();
    AnisotropicProperty(const QString& name, const QString& unit, const AnisotropicModel& model);
    AnisotropicProperty(const AnisotropicProperty&);

    AnisotropicModel& model();
    void setModel(const AnisotropicModel& model);

    QVector<DataType*> modelData() override;
    int setModelData(const QVector<DataType*>& data) override;

    void accept(MaterialPropertyVisitor* visitor) override;
};

class MaterialPropertyVisitor
{
public:
    virtual void visit(BasicProperty& p) = 0;
    virtual void visit(IsotropicProperty& p) = 0;
    virtual void visit(DiagonalProperty& p) = 0;
    virtual void visit(SymmetricProperty& p) = 0;
    virtual void visit(AnisotropicProperty& p) = 0;
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

    QString description() const;
    QString toString() const;
    void appendProperty(MaterialProperty* property);
    void removeProperty(int index);
    int modifyProperty(int index, MaterialProperty* property);
    MaterialProperty* property(int index) const;
};

#endif // MATERIAL_H
