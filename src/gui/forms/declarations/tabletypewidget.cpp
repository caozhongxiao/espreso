#include "tabletypewidget.h"

#include "../validators/validatordelegate.h"
#include "../validators/validatorfactory.h"

#include <QRegExp>
#include <QString>
#include <QStringList>
#include <QDebug>

QStringList TableTypeWidget::headlines()
{
    QStringList result;
    result << QObject::tr("x") << QObject::tr("f(x)");

    return result;
}

TableTypeWidget::TableTypeWidget(QWidget* parent) :
    TableWidget(2, TableTypeWidget::headlines(), parent)
{
    this->mTable->setItemDelegateForColumn(0, new ValidatorDelegate(new NumberValidatorFactory()));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new NumberValidatorFactory()));

    this->defaultValues << "" << "";
}

void TableTypeWidget::addData(const QString& data)
{
    QRegExp rxCase("case");
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
