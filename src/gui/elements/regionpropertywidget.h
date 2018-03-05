#ifndef REGIONPROPERTYWIDGET_H
#define REGIONPROPERTYWIDGET_H

#include <QWidget>
#include <QStandardItemModel>
#include <QTreeView>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/plane/planeelement.h"

#include "regionpairdialog.h"
#include "ecfobjecttreewidget.h"

namespace espreso
{

class RegionPropertyWidget : public ECFObjectTreeWidget
{
    Q_OBJECT

public:
    explicit RegionPropertyWidget(Mesh* mesh, PhysicsConfiguration* physics, const QString& label = "", QWidget *parent = 0);

    void addProperty(ECFObject* obj);

protected:
    virtual QDialog* createDialog(const QModelIndex& groupIndex, ECFParameter* param = nullptr) override;
    virtual QString dialogResult(QDialog* dialog) override;

    virtual void newItemAccepted(int, QString) override {}
    virtual void newItemRejected(int) override {}
    virtual void itemEditted(int, ECFParameter*) override {}

private:
    Mesh* m_mesh;
    PhysicsConfiguration* m_physics;
};

}

#endif // REGIONPROPERTYWIDGET_H
