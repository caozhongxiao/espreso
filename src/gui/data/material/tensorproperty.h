#ifndef TENSORPROPERTY_H
#define TENSORPROPERTY_H

#include "../namedentity.h"
#include "../datatype.h"
#include "model/isotropicmodel.h"
#include "model/diagonalmodel.h"
#include "model/symmetricmodel.h"
#include "model/anisotropicmodel.h"

class TensorProperty : public NamedEntity
{
private:
    QString mUnit;
    int mActiveModel;
    IsotropicModel mIM;
    DiagonalModel mDM;
    SymmetricModel mSM;
    AnisotropicModel mAM;

public:
    static const int MODEL_ISOTROPIC = 0x0001;
    static const int MODEL_DIAGONAL = 0x0002;
    static const int MODEL_SYMMETRIC = 0x0004;
    static const int MODEL_ANISOTROPIC = 0x0008;

    TensorProperty();
    TensorProperty(const QString& name, const QString& unit);
    TensorProperty(const TensorProperty& tp);

    int activeModel();
    void setActiveModel(int model);

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
