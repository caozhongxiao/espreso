#include "maptablewidget.h"

using namespace espreso;

MapTableWidget::MapTableWidget(const QStringList& headlines, QWidget *parent) :
    TableWidget(2, headlines, parent)
{

}

void MapTableWidget::addData(const QString& data)
{
    qWarning("MapTableWidget::addData: Empty method.");
}

QString MapTableWidget::data()
{
    return QLatin1String("");
}
