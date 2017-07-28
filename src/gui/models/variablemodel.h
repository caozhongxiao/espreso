#ifndef VARIABLEMODEL_H
#define VARIABLEMODEL_H

#include <QAbstractListModel>
#include <QVector>
#include "../data/datatype.h"
#include "../data/variable.h"

class VariableListModel : public QAbstractListModel
{
    Q_OBJECT
private:
    QVector<Variable> vars;

public:
    VariableListModel(QObject* parent = 0) : QAbstractListModel(parent) {}

    int rowCount(const QModelIndex &parent) const;
    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    Qt::ItemFlags flags(const QModelIndex &index) const;
    bool setData(const QModelIndex &index, const QVariant &value, int role);
    bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
    bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());
};

#endif // VARIABLEMODEL_H
