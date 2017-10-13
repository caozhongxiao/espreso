#ifndef DOUBLETABLEDELEGATE_H
#define DOUBLETABLEDELEGATE_H

#include <QItemDelegate>
#include "../data/common.h"

class DoubleTableDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    explicit DoubleTableDelegate(QObject* parent = 0);

protected:
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    void setEditorData(QWidget * editor, const QModelIndex & index) const override;
    void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
    void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;
};

#endif // DOUBLETABLEDELEGATE_H
