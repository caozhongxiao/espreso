#include "doubletabledelegate.h"
#include <QLineEdit>
#include <QIntValidator>

DoubleTableDelegate::DoubleTableDelegate(QObject *parent) : QItemDelegate(parent)
{
}

QWidget* DoubleTableDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    QLineEdit *editor = new QLineEdit(parent);
    QRegExpValidator* validator = new QRegExpValidator(QRegExp(REGEXPR_DOUBLE), parent);
    editor->setValidator(validator);
    return editor;
}


void DoubleTableDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value =index.model()->data(index, Qt::EditRole).toString();
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    line->setText(value);
}


void DoubleTableDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    QString value = line->text();
    model->setData(index, value);
}


void DoubleTableDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}
