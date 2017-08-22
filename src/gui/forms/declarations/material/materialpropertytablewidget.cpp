#include "materialpropertytablewidget.h"
#include "ui_materialpropertytablewidget.h"

#include <QComboBox>
#include <QLabel>
#include <QLineEdit>
#include "../../../data/common.h"
#include "../../validators/expressionvalidator.h"
#include "../../validators/validatorfactory.h"
#include "../tablewidget.h"
#include "../datatypeeditwidget.h"

MaterialPropertyTableWidget::MaterialPropertyTableWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MaterialPropertyTableWidget)
{
    ui->setupUi(this);
}

MaterialPropertyTableWidget::~MaterialPropertyTableWidget()
{
    delete ui;
}

void MaterialPropertyTableWidget::addRow(const QString& name, DataType* data,
                                         const QString& unit, const QString& abbrev)
{
    int row = ui->grid->rowCount();

    QLabel* lblName = new QLabel(name, this);
    QLabel* lblUnit = new QLabel(unit, this);
    QLabel* lblAbbrev = new QLabel(abbrev, this);

    DataTypeEditWidget* dataWidget = new DataTypeEditWidget(data, this);
    QComboBox* cmbBox = dataWidget->createComboBox(this);
    //TODO: Variables

    ui->grid->addWidget(lblName, row, 0);
    ui->grid->addWidget(cmbBox, row, 1);
    ui->grid->addWidget(dataWidget, row, 2);
    ui->grid->addWidget(lblUnit, row, 3);
    ui->grid->addWidget(lblAbbrev, row, 4);
}
