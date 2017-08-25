#ifndef SYMMETRICMODEL_H
#define SYMMETRICMODEL_H

#include "propertymodel.h"

class SymmetricModel : public PropertyModel
{

public:
    SymmetricModel();
    SymmetricModel(const QList<DataType*>& values);

    virtual void setMatrix(const QList<DataType*>& values) override;

private:
    void init();
};

#endif // SYMMETRICMODEL_H
