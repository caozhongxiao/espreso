#ifndef REGIONMATERIALSWIDGET_H
#define REGIONMATERIALSWIDGET_H

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/plane/planeelement.h"

#include <QWidget>

namespace espreso
{

namespace Ui {
class RegionMaterialsWidget;
}

class RegionMaterialsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit RegionMaterialsWidget(QWidget *parent = 0);
    ~RegionMaterialsWidget();

private:
    Ui::RegionMaterialsWidget *ui;
};

}

#endif // REGIONMATERIALSWIDGET_H
