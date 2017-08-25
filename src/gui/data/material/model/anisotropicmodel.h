#ifndef ANISOTROPICMODEL_H
#define ANISOTROPICMODEL_H

#include "propertymodel.h"

class AnisotropicModel : public PropertyModel
{

public:
    AnisotropicModel();
    AnisotropicModel(const QList<DataType*>& values);

    virtual void setMatrix(const QList<DataType*>& values) override;

private:
    void init();
};

#endif // ANISOTROPICMODEL_H
