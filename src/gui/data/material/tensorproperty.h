#ifndef TENSORPROPERTY_H
#define TENSORPROPERTY_H

#include "../namedentity.h"
#include "../datatype.h"

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

class TensorProperty : public NamedEntity
{
private:
    QString mUnit;
    IsotropicModel mIM;
    DiagonalModel mDM;
    SymmetricModel mSM;
    AnisotropicModel mAM;

public:
    TensorProperty();
    TensorProperty(const QString& name, const QString& unit);
    TensorProperty(const TensorProperty& tp);

    void setIsotropicModel(const IsotropicModel& im);
    void setDiagonalModel(const DiagonalModel& dm);
    void setSymmetricModel(const SymmetricModel& sm);
    void setAnisotropicModel(const AnisotropicModel& am);

    QString unit() const;
    IsotropicModel isotropicModel() const;
    DiagonalModel diagonalModel() const;
    SymmetricModel symmetricModel() const;
    AnisotropicModel anisotropicModel() const;
};

#endif // TENSORPROPERTY_H
