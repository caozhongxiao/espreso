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

namespace espreso
{

namespace Ui {
class RegionPropertyWidget;
}

class RegionPropertyWidget : public QWidget
{
    Q_OBJECT

public:
    explicit RegionPropertyWidget(Mesh* mesh, PhysicsConfiguration* physics, QWidget *parent = 0, const QString& label = "");
    ~RegionPropertyWidget();

    void addProperty(ECFObject* obj);

private slots:
    void onActionNew();
    void onActionEdit();
    void onActionDelete();

    void on_tree_customContextMenuRequested(const QPoint &pos);

private:
    Ui::RegionPropertyWidget *ui;

    Mesh* m_mesh;
    PhysicsConfiguration* m_physics;

    QStandardItemModel* m_model;
    QVector<QStandardItem*> m_groups;
    QVector<ECFObject*> m_objs;
    QStandardItem* m_root;

    QTreeView* m_view;

    QAction* m_action_new;
    QAction* m_action_edit;
    QAction* m_action_delete;

    QModelIndex selectedItem();
};

}

#endif // REGIONPROPERTYWIDGET_H
