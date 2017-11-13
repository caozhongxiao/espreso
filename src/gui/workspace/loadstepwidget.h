#ifndef LOADSTEPWIDGET_H
#define LOADSTEPWIDGET_H

#include "../elements/ecfobjectwidget.h"
#include "../elements/fieldhandler.h"
#include "../elements/regionpropertywidget.h"

namespace espreso
{

class LoadstepWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    explicit LoadstepWidget(size_t id, Mesh* mesh, ECFObject* physics, QWidget* parent = 0);

protected:
    QWidget* initContainer() override;
    void drawObject(ECFObject*) override;

private:
    ECFObject* m_physics;
    ECFParameter* m_loadstep;
    Mesh* m_mesh;

    RegionPropertyWidget* m_properties;

    QWidget* m_widget;
};

}

#endif // LOADSTEPWIDGET_H
