#include "regionpropertywidget.h"
#include "ui_regionpropertywidget.h"

#include <QMenu>
#include <QTreeView>

#include "regionpairdialog.h"

using namespace espreso;

RegionPropertyWidget::RegionPropertyWidget(Mesh* mesh, PhysicsConfiguration* physics, QWidget *parent, const QString& label) :
    QWidget(parent),
    ui(new Ui::RegionPropertyWidget)
{
    ui->setupUi(this);

    this->m_mesh = mesh;
    this->m_physics = physics;

    if (label == "")
        ui->label->hide();
    else
        ui->label->setText(label);

    this->m_model = new QStandardItemModel();
    this->m_root = this->m_model->invisibleRootItem();

    QTreeView* view = new QTreeView(this);
    view->setModel(m_model);
    view->setEditTriggers(QAbstractItemView::NoEditTriggers);
    view->setContextMenuPolicy(Qt::CustomContextMenu);
    view->setHeaderHidden(true);
    ui->verticalLayout->addWidget(view);
    connect(view, &QTreeView::customContextMenuRequested, this, &RegionPropertyWidget::on_tree_customContextMenuRequested);
    this->m_view = view;

    this->m_action_new = new QAction(tr("&New"), this);
    connect(this->m_action_new, &QAction::triggered, this, &RegionPropertyWidget::onActionNew);

    this->m_action_edit = new QAction(tr("&Edit"), this);
    connect(this->m_action_edit, &QAction::triggered, this, &RegionPropertyWidget::onActionEdit);

    this->m_action_delete = new QAction(tr("&Delete"), this);
    connect(this->m_action_delete, &QAction::triggered, this, &RegionPropertyWidget::onActionDelete);
}

RegionPropertyWidget::~RegionPropertyWidget()
{
    delete ui;
}

void RegionPropertyWidget::addProperty(ECFObject *obj)
{
    QStandardItem* group = new QStandardItem(QString::fromStdString(obj->name));
    this->m_root->appendRow(group);
    this->m_groups.append(group);
    this->m_objs.append(obj);

    for (auto it = obj->parameters.begin(); it != obj->parameters.end(); ++it)
    {
        group->appendRow(new QStandardItem(QString::fromStdString((*it)->name)));
    }
}

void RegionPropertyWidget::onActionNew()
{
    QModelIndex groupIndex = this->selectedItem();
    if ( !(groupIndex.isValid()) ) return;

    RegionPairDialog* dialog;

    if (m_objs[groupIndex.row()]->metadata.datatype.size() == 2)
    {
        if (m_objs[groupIndex.row()]->metadata.datatype[1] == ECFDataType::EXPRESSION)
            dialog = RegionPairDialog::createRegionExpression(m_objs[groupIndex.row()], m_mesh);

        if (m_objs[groupIndex.row()]->metadata.datatype[1] == ECFDataType::MATERIAL)
            dialog = RegionPairDialog::createRegionMaterial(m_objs[groupIndex.row()], m_mesh,
                    static_cast<ECFObject*>(m_physics->getParameter("materials")));
    }

    if (dialog->exec() == QDialog::Accepted)
    {
        QString region = dialog->region();
        QStandardItem* item = new QStandardItem(region);
        this->m_groups[groupIndex.row()]->appendRow(item);
    }
}

void RegionPropertyWidget::onActionEdit()
{

}

void RegionPropertyWidget::onActionDelete()
{

}

QModelIndex RegionPropertyWidget::selectedItem()
{
    QModelIndexList indexList = m_view->selectionModel()->selectedIndexes();
    if (!indexList.count())
        return QModelIndex();

    QModelIndex clicked = indexList.at(0);
    QModelIndex parent = clicked.parent();

    if (!parent.isValid())
    {
        parent = clicked;
    }

    return parent;
}

void RegionPropertyWidget::on_tree_customContextMenuRequested(const QPoint &pos)
{
    QMenu treeMenu(this);
    treeMenu.addAction(this->m_action_new);
    treeMenu.addAction(this->m_action_edit);
    treeMenu.addAction(this->m_action_delete);
    treeMenu.exec(m_view->mapToGlobal(pos));
}
