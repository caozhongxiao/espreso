#ifndef VARIABLE_H
#define VARIABLE_H

#include <QString>
#include "datatype.h"
#include <QMetaType>

class Variable
{
private:
    QString mName;
    DataType* mData;

public:
    Variable();
    Variable(const Variable& other);
    ~Variable();

    Variable(const QString& name, DataType* data);

    QString name() const;
    DataType* data() const;
    QString toString() const;
};

Q_DECLARE_METATYPE(Variable);

#endif // VARIABLE_H
