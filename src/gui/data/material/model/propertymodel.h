#ifndef PROPERTYMODEL_H
#define PROPERTYMODEL_H

#include <QVector>
#include <QList>
#include "../../datatype.h"

class PropertyModel
{
public:
    virtual ~PropertyModel();
    virtual void setMatrix(const QList<DataType*>& values) = 0;
    QList<DataType*> matrix() const;

protected:
    PropertyModel();
    PropertyModel(const QList<DataType*>& values);
    QList<DataType*> mMatrix;
    void reportError();
    void clearMatrix();
};

#endif // PROPERTYMODEL_H
