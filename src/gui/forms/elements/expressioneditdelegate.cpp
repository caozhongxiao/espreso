#include "expressioneditdelegate.h"

ExpressionEditDelegate::ExpressionEditDelegate(QObject *parent) : QItemDelegate(parent)
{
}

QWidget* ExpressionEditDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    ExpressionEdit *editor = new ExpressionEdit(parent);

    connect(editor, &ExpressionEdit::validStateChanged, this, &ExpressionEditDelegate::changeValidState);

    return editor;
}


void ExpressionEditDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    ExpressionEdit *line = static_cast<ExpressionEdit*>(editor);
    line->setText(value);
}


void ExpressionEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    ExpressionEdit *line = static_cast<ExpressionEdit*>(editor);
    QString value = line->text();
    model->setData(index, value);
}


void ExpressionEditDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}

void ExpressionEditDelegate::changeValidState(bool valid)
{
    emit validStateChanged(valid);
}
