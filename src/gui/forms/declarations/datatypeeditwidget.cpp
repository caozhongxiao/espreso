#include "datatypeeditwidget.h"
#include "ui_datatypeeditwidget.h"

#include "../validators/validatordelegate.h"
#include <QPair>
#include <QDebug>
#include <QRegularExpression>

DataTypeEditWidget::DataTypeEditWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DataTypeEditWidget)
{
    ui->setupUi(this);

    this->createUi();

    this->activeType = 0;
    this->uiExpression->show();
}

DataTypeEditWidget::DataTypeEditWidget(const ECFValue& data, QWidget *parent) :
    DataTypeEditWidget(parent)
{
    this->uiExpression->hide();

    QString content = QString::fromStdString(data.getValue());
    QRegularExpression tableregex("switch");
    QRegularExpression piecewiseregex("if");

    if (tableregex.match(content).hasMatch())
    {
        this->initTable(data);
        this->activeType = 1;
    }
    else if (piecewiseregex.match(content).hasMatch())
    {
        this->initPiecewise(data);
        this->activeType = 2;
    }
    else
    {
        this->initExpression(data);
        this->activeType = 0;
    }
}

DataTypeEditWidget::~DataTypeEditWidget()
{
    delete ui;
}

void DataTypeEditWidget::createUi()
{
    ExpressionEdit* function = new ExpressionEdit(this);
    function->hide();

    TableTypeWidget* table = new TableTypeWidget(this);
    table->hide();

    PiecewiseTypeWidget* piecewise = new PiecewiseTypeWidget(this);
    piecewise->hide();

    this->uiExpression = function;
    this->uiTable = table;
    this->uiPiecewise = piecewise;

    ui->layout->addWidget(uiExpression);
    ui->layout->addWidget(uiTable);
    ui->layout->addWidget(uiPiecewise);
}

void DataTypeEditWidget::initExpression(const ECFValue& et)
{
    this->uiExpression->setText(QString::fromStdString(et.getValue()));
    this->uiExpression->show();
}

void DataTypeEditWidget::initTable(const ECFValue& tt)
{
    this->uiTable->addData(QString::fromStdString(tt.getValue()));
    this->uiTable->show();
}

void DataTypeEditWidget::initPiecewise(const ECFValue& pft)
{
    this->uiPiecewise->addData(QString::fromStdString(pft.getValue()));
    this->uiPiecewise->show();
}

QComboBox* DataTypeEditWidget::createComboBox(QWidget *parent)
{
    QComboBox* box = new QComboBox(parent);
    box->addItems(DataTypeEditWidget::typeNames());
    box->setCurrentIndex(this->activeType);

    connect(box, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &DataTypeEditWidget::changeType);

    return box;
}

bool DataTypeEditWidget::isValid()
{
    switch (this->activeType)
    {
        case 0:
            return this->uiExpression->isValid();
        case 1:
            return this->uiTable->isValid();
        case 2:
            return this->uiPiecewise->isValid();
        default:
            qCritical() << "DataTypeEditWidget: Unknown DataType ID" << this->activeType;
            return false;
    }
}

void DataTypeEditWidget::changeType(int index)
{
    switch (index)
    {
        case 0:
            this->uiTable->hide();
            this->uiPiecewise->hide();
            this->uiExpression->show();
            break;
        case 1:
            this->uiExpression->hide();
            this->uiPiecewise->hide();
            this->uiTable->show();
            break;
        case 2:
            this->uiExpression->hide();
            this->uiTable->hide();
            this->uiPiecewise->show();
            break;
        default:
            qCritical() << "DataTypeEditWidget: Unknown DataType ID"
                      << index;
            return;
    }

    this->activeType = index;
}

QStringList DataTypeEditWidget::typeNames()
{
    QStringList result;
    result << QObject::tr("Expression")
           << QObject::tr("Table")
           << QObject::tr("Piecewise function");

    return result;
}
