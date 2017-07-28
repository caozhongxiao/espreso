#include "treeitem.h"

TreeItem::TreeItem(const QVector<QVariant> &data, TreeItem *parentItem)
{
    this->itemData = data;
    this->itemParent = parentItem;
}

TreeItem::~TreeItem()
{
    qDeleteAll(this->childItems);
}

void TreeItem::appendChild(TreeItem* child)
{
    this->childItems.append(child);
}

TreeItem* TreeItem::child(int row)
{
    return this->childItems.value(row);
}

int TreeItem::childCount() const
{
    return this->childItems.count();
}

int TreeItem::columnCount() const
{
    return this->itemData.count();
}

int TreeItem::row() const
{
    if (itemParent)
        return this->itemParent->childItems.indexOf(const_cast<TreeItem*>(this));

    return 0;
}

QVariant TreeItem::data(int column) const
{
    return this->itemData.value(column);
}

TreeItem* TreeItem::parentItem()
{
    return this->itemParent;
}

bool TreeItem::setData(int column, const QVariant &value)
{
    if (column < 0 || column >= itemData.size())
        return false;

    itemData[column] = value;
    return true;
}

bool TreeItem::insertChildren(int position, int count, int columns)
{
    if (position < 0 || position > childItems.size())
        return false;

    for (int row = 0; row < count; ++row) {
        QVector<QVariant> data(columns);
        TreeItem* item = new TreeItem(data, this);
        this->childItems.insert(position, item);
    }

    return true;
}

bool TreeItem::removeChildren(int position, int count)
{
    if (position < 0 || position + count > childItems.size())
        return false;

    for (int row = 0; row < count; ++row)
        delete childItems.takeAt(position);

    return true;
}

bool TreeItem::insertColumns(int position, int columns)
{
    if (position < 0 || position > itemData.size())
        return false;

    for (int column = 0; column < columns; ++column)
        this->itemData.insert(position, QVariant());

    foreach (TreeItem *child, childItems)
        child->insertColumns(position, columns);

    return true;
}

bool TreeItem::removeColumns(int position, int columns)
{
    if (position < 0 || position + columns > itemData.size())
        return false;

    for (int column = 0; column < columns; ++column)
        itemData.remove(position);

    foreach (TreeItem *child, childItems)
        child->removeColumns(position, columns);

    return true;
}

