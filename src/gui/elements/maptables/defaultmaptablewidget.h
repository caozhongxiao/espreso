#ifndef DEFAULTMAPTABLEWIDGET_H
#define DEFAULTMAPTABLEWIDGET_H

#include "../maptablewidget.h"

#include "../../../config/configuration.h"

namespace espreso
{

class DefaultMapTableWidget : public MapTableWidget
{
public:
    DefaultMapTableWidget(ECFObject* map, const QStringList& headlines, QWidget* parent = 0);

    QVector<QPair<QString, QString> > dataInRows() override;

    bool isValid() override;
    QString errorMessage() override;

    void save() override;

protected:
    ECFObject* m_map;

private:
    void setupColumn(int col, ECFDataType type);
};

}

#endif // DEFAULTMAPTABLEWIDGET_H
