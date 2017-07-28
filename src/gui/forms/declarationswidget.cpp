#include "declarationswidget.h"
#include "ui_declarationswidget.h"

#include "treemodel.h"
#include "variabledialog.h"
#include "variable.h"

DeclarationsWidget::DeclarationsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DeclarationsWidget)
{
    ui->setupUi(this);

    ui->DeclarationTree->setModel(new TreeModel(tr("Declarations")));
    this->setupTree();
}

DeclarationsWidget::~DeclarationsWidget()
{
    delete ui;
}

void DeclarationsWidget::setupTree()
{
    QAbstractItemModel* model = ui->DeclarationTree->model();
    model->insertRows(0, 3);
    QVariant vars(tr("Variables"));
    QVariant csystems(tr("Coordinate Systems"));
    QVariant mats(tr("Materials"));
    QModelIndex varIndex = model->index(0, 0);
    QModelIndex csIndex = model->index(1, 0);
    QModelIndex matsIndex = model->index(2, 0);
    model->setData(varIndex, vars);
    model->setData(csIndex, csystems);
    model->setData(matsIndex, mats);

}

void DeclarationsWidget::on_DeclarationTree_customContextMenuRequested(const QPoint &pos)
{
    QModelIndex index = ui->DeclarationTree->model()->index(pos.x(), pos.y());
    if (!index.isValid()) return;

}

void DeclarationsWidget::createActions()
{
    this->newItem = new QAction(tr("&New"), this);
    connect(this->newItem, &QAction::triggered, this, &DeclarationsWidget::treeNewItem);

    this->editItem = new QAction(tr("&Edit"), this);
    connect(this->editItem, &QAction::triggered, this, &DeclarationsWidget::treeEditItem);

    this->delItem = new QAction(tr("&Delete"), this);
    connect(this->delItem, &QAction::triggered, this, &DeclarationsWidget::treeDelItem);
}

void DeclarationsWidget::treeNewItem()
{
    QModelIndexList indexList = ui->DeclarationTree->selectionModel()->selectedIndexes();
    if (!indexList.count())
        return;
    QModelIndex index = indexList.at(0);
    QVariant data = index.data();
    if (data.canConvert<Variable>())
    {

    }
}

void DeclarationsWidget::treeEditItem()
{

}

void DeclarationsWidget::treeDelItem()
{

}
