#ifndef TENSORPROPERTY_H
#define TENSORPROPERTY_H

#include "../namedentity.h"
#include "../datatype.h"
#include "tensorpropertymodel.h"

#include <QStringList>
#include <QObject>
#include <QVector>

class TensorProperty : public NamedEntity
{

public:
    TensorProperty();
    TensorProperty(const QString& name, const QString& unit);
    TensorProperty(const TensorProperty& tp);

    int activeModel();
    void setActiveModel(int model);

    void appendModel(const TensorPropertyModel& model);
    const TensorPropertyModel& model(int index) const;

    auto modelBegin();
    auto modelEnd();

    QString unit() const;\

private:
    QString mUnit;
    int mActiveModel;
    QVector<TensorPropertyModel> mModels;

};

#endif // TENSORPROPERTY_H
