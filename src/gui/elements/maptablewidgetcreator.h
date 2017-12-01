#ifndef MAPTABLEWIDGETCREATOR_H
#define MAPTABLEWIDGETCREATOR_H

#include "maptablewidget.h"
#include "maptables/defaultmaptablewidget.h"
#include "maptables/boolmaptablewidget.h"

namespace espreso
{

class MapTableWidgetCreator
{
public:
    MapTableWidgetCreator();
    MapTableWidget* create(ECFObject* map, QWidget* parent = 0);

private:
    QVector<ECFDataType> m_first;
    QVector<ECFDataType> m_second;
};

}

#endif // MAPTABLEWIDGETCREATOR_H
