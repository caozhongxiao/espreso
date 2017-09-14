#include "expressioneditdelegate.h"

using namespace espreso;

ExpressionEditDelegate::ExpressionEditDelegate(QObject *parent) : QItemDelegate(parent)
{
}

QWidget* ExpressionEditDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    QLineEdit *editor = new QLineEdit(parent);

    return editor;
}


void ExpressionEditDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    line->setText(value);
}


void ExpressionEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    QString value = line->text();
    model->setData(index, value);
}


void ExpressionEditDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}

void ExpressionEditDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    QVariant data = index.data();
    QString text = data.toString();
    QStyleOptionViewItem newOption(option);
    if (!ExpressionEdit::validate(text))
    {
        newOption.font.setBold(true);
        newOption.palette.setColor(QPalette::Text, Qt::red);
    }

    QItemDelegate::paint(painter, newOption, index);
}
