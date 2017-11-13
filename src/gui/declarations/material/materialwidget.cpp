#include "materialwidget.h"

#include <QVBoxLayout>
#include <QScrollArea>
#include <QComboBox>
#include <QMessageBox>
#include <QLabel>
#include <QScrollArea>
#include <QLineEdit>
#include <QDebug>
#include "../../elements/optionhandler.h"
#include "../../elements/formwidget.h"
#include "../../elements/boolhandler.h"
#include "materialpropertytablewidget.h"

using namespace espreso;

MaterialWidget::MaterialWidget(MaterialConfiguration* material,
                               const QVector<std::string>& materialNames,
                               QWidget *parent) :
    ECFObjectWidget(material, parent)
{
    this->m_names = materialNames;
}


QWidget* MaterialWidget::initContainer()
{
    QScrollArea* area = new QScrollArea;

    QWidget* widget = new QWidget(area);
    this->m_widget = widget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    widget->setLayout(w_layout);

    this->createHeadline(this->m_obj, widget);

    area->setWidgetResizable(true);
    area->setWidget(widget);

    return area;
}

void MaterialWidget::drawObject(ECFObject* obj)
{
    QWidget* widget = new QWidget;
    QVBoxLayout* layout = new QVBoxLayout;
    widget->setLayout(layout);

    this->m_widget->layout()->addWidget(widget);

    this->createHeadline(obj, widget);

    this->processParameters(obj, widget);

    QSpacerItem* verticalSpacer = new QSpacerItem(0,
                                                  0,
                                                  QSizePolicy::Minimum,
                                                  QSizePolicy::Expanding);
    layout->addItem(verticalSpacer);
}

bool MaterialWidget::isValid()
{
    if (!ECFObjectWidget::isValid())
        return false;

    ISavableObject* info = m_savables.first();
    info->saveState();
    info->save();

    std::string newName =
            QString::fromStdString(
                m_obj->getParameter(std::string("name"))->getValue()
                )
                .toUpper()
                .toStdString();

    if (m_names.indexOf(newName) != -1)
    {
        QMessageBox msg;
        msg.setWindowTitle(tr("Error"));
        msg.setText(tr("Material with same name already exists!"));
        msg.exec();

        info->restoreState();

        return false;
    }

    foreach (ISavableObject* obj, m_savables) {
        obj->save();
    }

    return true;
}
