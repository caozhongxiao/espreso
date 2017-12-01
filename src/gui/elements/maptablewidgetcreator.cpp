#include "maptablewidgetcreator.h"

using namespace espreso;

MapTableWidgetCreator::MapTableWidgetCreator()
{
    this->m_first << ECFDataType::POSITIVE_INTEGER
            << ECFDataType::NONNEGATIVE_INTEGER
            << ECFDataType::STRING;
    this->m_second << ECFDataType::POSITIVE_INTEGER
             << ECFDataType::NONNEGATIVE_INTEGER
             << ECFDataType::INTERVAL;
}

MapTableWidget* MapTableWidgetCreator::create(ECFObject* map, QWidget* parent)
{
    if (map->metadata.datatype.size() != 2)
    {
        qCritical("MapTableWidget: Input map should contain two datatypes!");
        return nullptr;
    }
    if (map->metadata.description.size() != 2)
    {
        qCritical("MapTableWidget: Input map should contain two descriptions!");
        return nullptr;
    }

    QStringList headlines;
    headlines << QString::fromStdString(map->metadata.description[0])
            << QString::fromStdString(map->metadata.description[1]);

    if (this->m_first.contains(map->metadata.datatype[0])
            && this->m_second.contains(map->metadata.datatype[1]))
    {
        return new DefaultMapTableWidget(map, headlines, parent);
    }
    else if (this->m_first.contains(map->metadata.datatype[0])
             && map->metadata.datatype[1] == ECFDataType::BOOL)
    {
        return new BoolMapTableWidget(map, headlines, parent);
    }

    return nullptr;
}
