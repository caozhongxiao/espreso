#include "piecewisetypewidget.h"

#include "../elements/expressioneditdelegate.h"
#include "../elements/expressionedit.h"
#include "../validators/validatordelegate.h"
#include "../elements/comboboxdelegate.h"

#include <QDebug>

QStringList PiecewiseTypeWidget::headlines()
{
    QStringList result;
    result << QString("")
           << QObject::tr("Lower bound")
           << QObject::tr("Upper bound")
           << QString("")
           << QObject::tr("f(x)");

    return result;
}

PiecewiseTypeWidget::PiecewiseTypeWidget(QWidget* parent) :
    TableWidget(5, PiecewiseTypeWidget::headlines(), parent)
{
    QStringList leftCmbOptions;
    leftCmbOptions << "<" << "(";
    QStringList rightCmbOptions;
    rightCmbOptions << ">" << ")";
    this->mTable->setItemDelegateForColumn(0, new ComboBoxDelegate(leftCmbOptions, this));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new NumberValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(2, new ValidatorDelegate(new NumberValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(3, new ComboBoxDelegate(rightCmbOptions, this));
    this->mTable->setItemDelegateForColumn(4, new ExpressionEditDelegate(this));

    this->defaultValues << "<" << "" << "" << ")" << "";
}

bool PiecewiseTypeWidget::isValid()
{
    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QModelIndex index = this->mModel->index(row, 4);
        QString expression = index.data().toString();
        if (!ExpressionEdit::validate(expression))
            return false;
    }

    return TableWidget::isValid();
}


void PiecewiseTypeWidget::addData(const QString& data)
{
    QRegExp rxCase("if");
    QStringList ifs = data.split(rxCase, QString::SkipEmptyParts);

    foreach (QString v_if, ifs) {

        QVector<QString> row;

        QStringList middle = v_if.split(')', QString::SkipEmptyParts);

        QString condition = middle.at(0);
        condition.remove(0, 1);
        QStringList bounds = condition.split("and", QString::SkipEmptyParts);

        QString leftBound = bounds.at(0);
        QString rightBound = bounds.at(1);

        qDebug() << leftBound;
        qDebug() << rightBound;

        QString left;
        if (leftBound.contains(">="))
        {
            left = leftBound.split(">=", QString::SkipEmptyParts).at(1);
            row.append("<");
        }
        else if (leftBound.contains('>'))
        {
            left = leftBound.split(">", QString::SkipEmptyParts).at(1);
            row.append("(");
        }
        row.append(left.trimmed());

        QString right;
        if (rightBound.contains("<="))
        {
            right = rightBound.split("<=", QString::SkipEmptyParts).at(1);
            row.append(right.trimmed());
            row.append(">");
        }
        else if (rightBound.contains('<'))
        {
            right = rightBound.split('<', QString::SkipEmptyParts).at(1);
            row.append(right.trimmed());
            row.append(")");
        }

        QString expr = middle.at(1);
        expr = expr.split(';', QString::SkipEmptyParts).at(0);
        expr = expr.trimmed();
        row.append(expr);

        this->addRow(row);
    }
}

QString PiecewiseTypeWidget::columnDefaultValue(int column) const
{
    return this->defaultValues.at(column);
}
