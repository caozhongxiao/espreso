#include "materialpropertytablewidget.h"
#include "ui_materialpropertytablewidget.h"

#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

#include "../../../data/common.h"
#include "../../validators/validatorfactory.h"
#include "../tablewidget.h"
#include "../datatypeeditwidget.h"

using namespace espreso;

MaterialPropertyTableWidget::MaterialPropertyTableWidget(QWidget *parent, bool withHeader) :
    QWidget(parent),
    ui(new Ui::MaterialPropertyTableWidget)
{
    ui->setupUi(this);

    if (withHeader) this->createHeader();
}

MaterialPropertyTableWidget::~MaterialPropertyTableWidget()
{
    delete ui;
}

void MaterialPropertyTableWidget::createHeader()
{
    QLabel* lblCol1 = new QLabel(tr("Name"));
    QLabel* lblCol2 = new QLabel(tr("Type"));
    QLabel* lblCol3 = new QLabel(tr("Definition"));
    QLabel* lblCol4 = new QLabel(tr("Unit"));
    QLabel* lblCol5 = new QLabel(tr("Symbol"));

    ui->grid->addWidget(lblCol1, 0, 0);
    ui->grid->addWidget(lblCol2, 0, 1);
    ui->grid->addWidget(lblCol3, 0, 2);
    ui->grid->addWidget(lblCol4, 0, 3);
    ui->grid->addWidget(lblCol5, 0, 4);
}

void MaterialPropertyTableWidget::addProperty(const ECFParameter& property)
{
    this->addRow(QString::fromStdString(property.metadata.description.at(0)),
                 property,
                 QString::fromStdString(property.metadata.unit),
                 QString::fromStdString(property.name));
}

void MaterialPropertyTableWidget::addRow(const QString& name, const ECFParameter& data,
                                         const QString& unit, const QString& symbol)
{
    int row = ui->grid->rowCount();

    QLabel* lblName = new QLabel(name, this);
    QLabel* lblUnit = new QLabel(unit, this);
    QLabel* lblSymbol = new QLabel(symbol, this);

    DataTypeEditWidget* dataWidget = new DataTypeEditWidget(data, this);
    QComboBox* cmbBox = dataWidget->createComboBox(this);
    //TODO: Variables

    Qt::Alignment alignment = Qt::AlignHCenter | Qt::AlignTop;
    ui->grid->addWidget(lblName, row, 0, alignment);
    ui->grid->addWidget(cmbBox, row, 1, alignment);
    ui->grid->addWidget(dataWidget, row, 2);
    ui->grid->addWidget(lblUnit, row, 3, alignment);
    ui->grid->addWidget(lblSymbol, row, 4, alignment);
}
