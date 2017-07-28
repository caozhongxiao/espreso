#include "variablemodel.h"

int VariableListModel::rowCount(const QModelIndex &parent) const
{
    return this->vars.size();
}

QVariant VariableListModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (index.row() >= this->rowCount(index))
        return QVariant();

    if (role == Qt::DisplayRole)
    {
        return this->vars[index.row()].toString();
    }
    else if (role ==  Qt::EditRole)
    {
        Variable var = this->vars[index.row()];
        QVariant result;
        result.setValue(var);
        return result;
    }
    else
        return QVariant();
}

QVariant VariableListModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role != Qt::DisplayRole)
        return QVariant();

    if (orientation == Qt::Horizontal)
        return QString(tr("Column %1")).arg(section);
    else
        return QString(tr("Row %1")).arg(section);
}

Qt::ItemFlags VariableListModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return Qt::ItemIsEnabled;

    return QAbstractItemModel::flags(index);
}

bool VariableListModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (index.isValid() && role == Qt::EditRole) {
        this->vars.replace(index.row(), value.value<Variable>());
        emit dataChanged(index, index);
        return true;
    }
    return false;
}

bool VariableListModel::insertRows(int position, int rows, const QModelIndex &index)
{
    beginInsertRows(QModelIndex(), position, position+rows-1);

    for (int row = 0; row < rows; ++row) {
        Variable v;
        this->vars.append(v);
    }

    endInsertRows();
    return true;
}

bool VariableListModel::removeRows(int position, int rows, const QModelIndex &index)
{
    beginRemoveRows(QModelIndex(), position, position+rows-1);

    for (int row = 0; row < rows; ++row) {
        this->vars.removeAt(position);
    }

    endRemoveRows();
    return true;
}
