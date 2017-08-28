#ifndef TENSORPROPERTYMODEL_H
#define TENSORPROPERTYMODEL_H

#include <QVector>
#include <QList>
#include <QString>
#include "materialproperty.h"

class TensorPropertyModelItem : public MaterialProperty
{
public:
    TensorPropertyModelItem() : MaterialProperty() {}
    TensorPropertyModelItem(const QString& name, const QString& unit, const QString& symbol,
                   DataType* data) :
        MaterialProperty(name, unit, symbol, data) {}
    TensorPropertyModelItem(const TensorPropertyModelItem& sp) : MaterialProperty(sp) {}
};


class TensorPropertyModel : public NamedEntity
{
public:
    TensorPropertyModel();
    TensorPropertyModel(int size, const QString& name, const QList<TensorPropertyModelItem>& values);
    TensorPropertyModel(const TensorPropertyModel& other);
    ~TensorPropertyModel();

    const QVector<TensorPropertyModelItem>& items() const;

    void setCellValue(int row, int column, TensorPropertyModelItem value);
    const TensorPropertyModelItem& cellValue(int row, int column) const;

    int size() const;

private:
    QVector<TensorPropertyModelItem> mMatrix;
    int mSize;
};

#endif // TENSORPROPERTYMODEL_H
