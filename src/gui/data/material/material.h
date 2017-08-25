#ifndef MATERIAL_H
#define MATERIAL_H

#include <QString>
#include "../namedentity.h"
#include "scalarproperty.h"
#include "tensorproperty.h"

class Material : public NamedEntity
{
private:
    QString mDescription;
    QVector<ScalarProperty> mScalars;
    QVector<TensorProperty> mTensors;

public:
    Material();
    Material(const QString& name, const QString& description);
    Material(const Material&);

    QString description() const;
    int append(const ScalarProperty& property);
    int append(const TensorProperty& property);
    void removeScalar(int index);
    void removeTensor(int index);
    const ScalarProperty& scalar(int index) const;
    const TensorProperty& tensor(int index) const;
};

#endif // MATERIAL_H
