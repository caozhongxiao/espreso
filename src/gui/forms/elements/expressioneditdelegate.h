#ifndef EXPRESSIONEDITDELEGATE_H
#define EXPRESSIONEDITDELEGATE_H

#include <QItemDelegate>
#include "expressionedit.h"

class ExpressionEditDelegate: public QItemDelegate
{
    Q_OBJECT

public:
    ExpressionEditDelegate(QObject* parent = nullptr);

signals:
    void validStateChanged(bool valid);

private slots:
    void changeValidState(bool valid);

protected:
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    void setEditorData(QWidget * editor, const QModelIndex & index) const override;
    void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
    void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;
};

#endif // EXPRESSIONEDITDELEGATE_H
