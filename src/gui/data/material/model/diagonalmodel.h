#ifndef DIAGONALMODEL_H
#define DIAGONALMODEL_H

#include "propertymodel.h"

class DiagonalModel : public PropertyModel
{

public:
    DiagonalModel();
    DiagonalModel(const QList<DataType*>& values);

    virtual void setMatrix(const QList<DataType*>& values) override;

private:
    void init();
};

#endif // DIAGONALMODEL_H
