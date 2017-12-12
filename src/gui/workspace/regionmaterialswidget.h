#ifndef REGIONMATERIALSWIDGET_H
#define REGIONMATERIALSWIDGET_H

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/plane/planeelement.h"

#include "../elements/isavableobject.h"
#include "../elements/ivalidatableobject.h"

#include <QWidget>

namespace espreso
{

namespace Ui {
class RegionMaterialsWidget;
}

class RegionMaterialsWidget : public QWidget, public ISavableObject, public IValidatableObject
{
    Q_OBJECT

    virtual void save() override {}
    virtual bool isValid() override { return true; }
    virtual QString errorMessage() override { return QLatin1String(""); }

public:
    explicit RegionMaterialsWidget(Mesh* mesh, PhysicsConfiguration* physics, QWidget *parent = 0);
    ~RegionMaterialsWidget();

private:
    Ui::RegionMaterialsWidget *ui;
};

}

#endif // REGIONMATERIALSWIDGET_H
