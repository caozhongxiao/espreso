#ifndef MAPTABLEWIDGET_H
#define MAPTABLEWIDGET_H

#include "../declarations/tablewidget.h"
#include "isavableobject.h"

namespace espreso
{

class MapTableWidget : public TableWidget, public ISavableObject
{
    Q_OBJECT

public:
    virtual QVector<QPair<QString, QString> > dataInRows() = 0;

    virtual void addData(const QString& data) override;
    virtual QString data() override;

protected:
    MapTableWidget(const QStringList& headlines, QWidget *parent = 0);
};

}

#endif // MAPTABLEWIDGET_H
