#ifndef PIECEWISETYPEWIDGET_H
#define PIECEWISETYPEWIDGET_H

#include "tablewidget.h"
#include "../../data/datatype.h"

class PiecewiseTypeWidget : public TableWidget
{
    Q_OBJECT

public:
    static QStringList headlines();

    PiecewiseTypeWidget(QWidget* parent = 0);
    bool isValid() override;

    virtual void addData(const QString&) override;

protected:
    virtual QString columnDefaultValue(int column) const override;

private:
    QVector<QString> defaultValues;
};

#endif // PIECEWISETYPEWIDGET_H
