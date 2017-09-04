#ifndef TABLETYPEWIDGET_H
#define TABLETYPEWIDGET_H

#include "tablewidget.h"
#include "../../data/datatype.h"

class TableTypeWidget : public TableWidget
{
    Q_OBJECT

public:
    static QStringList headlines();

    TableTypeWidget(QWidget *parent = 0);

    virtual void addData(const QString& data) override;

protected:
    virtual QString columnDefaultValue(int column) const override;

private:
    QVector<QString> defaultValues;
};

#endif // TABLETYPEWIDGET_H
