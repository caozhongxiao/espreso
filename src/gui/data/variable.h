#ifndef VARIABLE_H
#define VARIABLE_H


class Variable
{
private:
    QString name;
    DataType* data;

public:
    Variable();
    Variable(const Variable &other);
    ~Variable();

    Variable(QString name, DataType* data);

    QString name() const;
    const DataType* data() const;
    QString toString() const;
};
Q_DECLARE_METATYPE(Variable);
#endif // VARIABLE_H
