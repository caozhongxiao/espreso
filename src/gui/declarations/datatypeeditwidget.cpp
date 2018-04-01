#include "datatypeeditwidget.h"
#include "ui_datatypeeditwidget.h"

#include "../validators/validatordelegate.h"
#include <QPair>
#include <QDebug>
#include <QRegularExpression>

using namespace espreso;

DataTypeEditWidget::DataTypeEditWidget(QWidget *parent) :
    TextItemWidget(parent),
    ui(new Ui::DataTypeEditWidget)
{
    ui->setupUi(this);

    this->m_cmb = this->createComboBox();
    ui->layout->addWidget(m_cmb);
    this->m_cmb->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    this->m_cmb->hide();
}

DataTypeEditWidget::DataTypeEditWidget(const std::vector<std::string>& variables, QWidget* parent) :
    DataTypeEditWidget(parent)
{
    this->m_param = nullptr;
    this->m_param_variables = variables;

    this->createUi();
    this->initExpression();
    this->activeType = 0;
}

DataTypeEditWidget::DataTypeEditWidget(ECFParameter* data, QWidget *parent) :
    DataTypeEditWidget(parent)
{    
    this->m_param = data;
    this->createUi();

    QString content = QString::fromStdString(data->getValue());
	this->parseValue(content);
}

DataTypeEditWidget::~DataTypeEditWidget()
{
    delete ui;
}

void DataTypeEditWidget::createUi()
{
    ExpressionEdit* function = new ExpressionEdit(
                "",
                param_variables(),
                this);
    function->hide();

	TableTypeWidget* table = new TableTypeWidget(this);
    table->hide();

    PiecewiseTypeWidget* piecewise = new PiecewiseTypeWidget(
                param_variables(),
                this);
    piecewise->hide();

    this->uiExpression = function;
    this->uiTable = table;
    this->uiPiecewise = piecewise;

    ui->layout->addWidget(uiExpression);
    ui->layout->addWidget(uiTable);
    ui->layout->addWidget(uiPiecewise);
}

void DataTypeEditWidget::parseValue(const QString &value)
{
	QRegularExpression tableregex("tabular|TABULAR");
	QRegularExpression piecewiseregex("if|IF");

	if (tableregex.match(value).hasMatch())
	{
		this->initTable();
		this->activeType = 1;
	}
	else if (piecewiseregex.match(value).hasMatch())
	{
		this->initPiecewise();
		this->activeType = 2;
	}
	else
	{
		this->initExpression();
		this->activeType = 0;
	}
}

void DataTypeEditWidget::initExpression()
{
    this->m_cmb->setCurrentIndex(0);
    this->uiExpression->setText(param_getValue());
    this->uiExpression->show();
}

void DataTypeEditWidget::initTable()
{
    this->m_cmb->setCurrentIndex(1);
    this->uiTable->addData(param_getValue());
    this->uiTable->show();
}

void DataTypeEditWidget::initPiecewise()
{
    this->m_cmb->setCurrentIndex(2);
    this->uiPiecewise->addData(param_getValue());
    this->uiPiecewise->show();
}

QString DataTypeEditWidget::param_getValue()
{
	if (!this->m_text.isEmpty()) return this->m_text;
    if (this->m_param != nullptr) return QString::fromStdString(m_param->getValue());
    return QString();
}

std::vector<std::string> DataTypeEditWidget::param_variables()
{
    if (this->m_param != nullptr) return m_param->metadata.variables;
    return this->m_param_variables;
}

QComboBox* DataTypeEditWidget::createComboBox(QWidget *parent)
{
    QComboBox* box = new QComboBox(parent);
    box->addItems(DataTypeEditWidget::typeNames());
    box->setCurrentIndex(this->activeType);

    connect(box, SIGNAL(currentIndexChanged(int)),
            this, SLOT(changeType(int)));

    return box;
}

void DataTypeEditWidget::setComboBox(bool show)
{
    if (show)
        this->m_cmb->show();
    else
        this->m_cmb->hide();
}

void DataTypeEditWidget::setSharedData(DataTypeEditWidgetFactoryData *data)
{
	this->m_shared = data;
	this->m_cmb->setCurrentIndex(data->type);
	this->changeType(data->type);
}

int DataTypeEditWidget::datatype()
{
    return this->activeType;
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

QString DataTypeEditWidget::errorMessage()
{
    switch (this->activeType)
    {
        case 0:
            return this->uiExpression->errorMessage();
        case 1:
            return this->uiTable->errorMessage();
        case 2:
            return this->uiPiecewise->errorMessage();
        default:
            qCritical() << "DataTypeEditWidget: Unknown DataType ID" << this->activeType;
            return QLatin1String("");
    }
}

void DataTypeEditWidget::save()
{
    if (this->m_param == nullptr) return;

    if (this->activeType == 0)
    {
        this->m_param->setValue(uiExpression->text().toStdString());
    }
    else if (this->activeType == 1)
    {
        this->m_param->setValue(uiTable->data().toStdString());
    }
    else if (this->activeType == 2)
    {
        this->m_param->setValue(uiPiecewise->data().toStdString());
    }
}

void DataTypeEditWidget::setText(const QString &text)
{
	this->m_text = text;
	if (!text.isEmpty())
	{
		this->parseValue(text);
		if (this->m_shared == nullptr) return;
		this->m_shared->valid = this->isValid();
		this->m_shared->error_message = this->errorMessage();
	}
}

QString DataTypeEditWidget::text()
{
	return this->value();
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
	if (this->m_shared != nullptr) this->m_shared->type = index;
    emit finished(this);
}

QStringList DataTypeEditWidget::typeNames()
{
    QStringList result;
    result << QObject::tr("Expression")
           << QObject::tr("Table")
           << QObject::tr("Piecewise function");

    return result;
}

QString DataTypeEditWidget::value()
{
    if (this->activeType == 0)
    {
        return uiExpression->text();
    }
    else if (this->activeType == 1)
    {
        return uiTable->data();
    }
    else if (this->activeType == 2)
    {
        return uiPiecewise->data();
    }
    else
    {
        return "";
    }
}

void DataTypeEditWidget::setValue(const QString &value)
{
    QRegularExpression tableregex("switch|SWITCH");
    QRegularExpression piecewiseregex("if|IF");

    if (tableregex.match(value).hasMatch())
    {
        this->uiTable->addData(value);
        this->changeType(1);
    }
    else if (piecewiseregex.match(value).hasMatch())
    {
        this->uiPiecewise->addData(value);
        this->changeType(2);
    }
    else
    {
        this->uiExpression->setText(value);
        this->changeType(0);
    }
}
