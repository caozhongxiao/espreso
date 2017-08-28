#ifndef TENSORPROPERTYMODEL_H
#define TENSORPROPERTYMODEL_H

#include <QVector>
#include <QList>
#include <QString>
#include "../datatype.h"
#include "../namedentity.h"

class TensorPropertyModel : public NamedEntity
{
public:
    TensorPropertyModel();
    TensorPropertyModel(int size, const QString& name);
    TensorPropertyModel(int size, const QString& name, const QList<DataType*>& values);
    TensorPropertyModel(const TensorPropertyModel& other);
    ~TensorPropertyModel();

    void setMatrix(const QList<DataType*>& values);
    QVector<DataType*> matrix() const;

    void setCellValue(int row, int column, DataType* value);
    DataType* cellValue(int row, int column) const;

    int size() const;

private:
    QVector<DataType*> mMatrix;
    int mSize;

    void clearMatrix();
};

#endif // TENSORPROPERTYMODEL_H
