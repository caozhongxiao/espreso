#ifndef ISOTROPICMODEL_H
#define ISOTROPICMODEL_H

#include "propertymodel.h"

class IsotropicModel : public PropertyModel
{

public:
    IsotropicModel();
    IsotropicModel(const QList<DataType*>& values);

    virtual void setMatrix(const QList<DataType*>& values) override;

private:
    void init();
};

#endif // ISOTROPICMODEL_H
