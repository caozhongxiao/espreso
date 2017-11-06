#include "ecfobjectwidget.h"
#include "ui_ecfobjectwidget.h"

#include <QMessageBox>
#include <QScrollArea>

using namespace espreso;

ECFObjectWidget::ECFObjectWidget(ECFObject* object, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ECFObjectWidget)
{
    this->m_obj = object;
    ui->setupUi(this);
}

ECFObjectWidget::~ECFObjectWidget()
{
    delete ui;
}

void ECFObjectWidget::init()
{
    this->drawMe();
}

void ECFObjectWidget::drawMe()
{
    this->m_savables.clear();
    this->m_validatables.clear();

    this->m_container = this->initContainer();
    ui->layout->addWidget(this->m_container);

    this->drawObject(this->m_obj);
}

OptionHandler* ECFObjectWidget::createOption(ECFParameter* option, QWidget* parent, bool withLabel)
{
    OptionHandler* handler = new OptionHandler(option, parent, withLabel);
    connect(handler, &OptionHandler::optionChanged, this, &ECFObjectWidget::redraw);

    return handler;
}

BoolHandler* ECFObjectWidget::createBool(ECFParameter* param, QWidget* parent)
{
    BoolHandler* handler = new BoolHandler(param, parent);
    connect(handler, &BoolHandler::stateChanged, this, &ECFObjectWidget::redraw);

    return handler;
}

void ECFObjectWidget::redraw()
{
    foreach (ISavableObject* obj, m_savables) {
        obj->save();
    }

    ui->layout->removeWidget(this->m_container);

    this->drawMe();
}

bool ECFObjectWidget::validate()
{
    foreach (IValidatableObject* obj, m_validatables) {
        if (!obj->isValid())
        {
            QMessageBox msg;
            msg.setWindowTitle(tr("Error"));
            msg.setText(obj->errorMessage());
            msg.exec();
            return false;
        }
    }

    return true;
}

bool ECFObjectWidget::isValid()
{
    foreach (IValidatableObject* obj, m_validatables) {
        if (!obj->isValid())
        {
            this->m_errormsg = obj->errorMessage();
            return false;
        }
    }

    return true;
}

QString ECFObjectWidget::errorMessage()
{
    this->isValid();

    return this->m_errormsg;
}

void ECFObjectWidget::save()
{
    foreach (ISavableObject* obj, m_savables)
    {
        obj->save();
    }
}
