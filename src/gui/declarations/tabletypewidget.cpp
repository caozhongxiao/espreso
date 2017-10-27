#include "tabletypewidget.h"

#include "../validators/validatordelegate.h"
#include "../validators/validatorfactory.h"

#include <QRegExp>
#include <QString>
#include <QStringList>
#include <QDebug>

using namespace espreso;

QStringList TableTypeWidget::headlines()
{
    QStringList result;
    result << QObject::tr("x") << QObject::tr("f(x)");

    return result;
}

TableTypeWidget::TableTypeWidget(QWidget* parent) :
    TableWidget(2, TableTypeWidget::headlines(), parent)
{
    this->mTable->setItemDelegateForColumn(0, new ValidatorDelegate(new DoubleValidatorFactory()));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new DoubleValidatorFactory()));

    this->defaultValues << "" << "";
}

void TableTypeWidget::addData(const QString& data)
{
    QRegExp rxCase("CASE");
    QStringList cases = data.split(rxCase, QString::SkipEmptyParts);
    cases.removeFirst();

    foreach (QString v_case, cases) {

            QStringList middle =  v_case.split(':', QString::SkipEmptyParts);

            QString condition = middle.at(0);
            QString middle_fx = middle.at(1);
            QString m_fx = middle_fx.split(';', QString::SkipEmptyParts).at(0);

            QStringList condition_parts = condition.split(QRegExp("=="), QString::SkipEmptyParts);
            QString m_x = condition_parts.at(1);

            QString x = m_x.simplified();
            x.replace(" ", "");

            QString fx = m_fx.simplified();
            fx.replace(" ", "");

            QVector<QString> pair;
            pair << x << fx;
            this->addRow(pair);
    }
}

QString TableTypeWidget::columnDefaultValue(int column) const
{
    return this->defaultValues.at(column);
}

QString TableTypeWidget::data()
{
    QLatin1String start("SWITCH {");
    QLatin1String end("DEFAULT: 0;}");

    QStringList cases;
    for (int i = 0; i < mModel->rowCount() - 1; ++i)
    {
        QLatin1String _case("CASE X == ");
        QModelIndex indexX = mModel->index(i, 0);
        QString _x = mModel->data(indexX).toString();
        QModelIndex indexY = mModel->index(i, 1);
        QString _y = mModel->data(indexY).toString();
        QString line = _case + _x + QLatin1String(" : ") + _y + QLatin1String("; ");
        cases << line;
    }

    return start + cases.join(QLatin1String(" ")) + end;
}
