#include "integertabwidget.h"

#include <QLabel>
#include <QPushButton>

using namespace espreso;

IntegerTabWidget::IntegerTabWidget(ECFObject* map,
                                   std::unique_ptr<ECFObjectWidgetFactory> factory,
                                   QWidget* parent) :
    QWidget(parent)
{
    this->m_map = map;
    this->m_factory = std::move(factory);

    QVBoxLayout* layout = new QVBoxLayout();
    this->setLayout(layout);

    QHBoxLayout* header = new QHBoxLayout();

    QLabel* lbl = new QLabel(QString::fromStdString(this->m_map->metadata.description.at(0)));
    header->addWidget(lbl);

    QPushButton* btnAdd = new QPushButton(tr("Add"));
    header->addWidget(btnAdd);
    connect(btnAdd, &QPushButton::pressed,
            this, &IntegerTabWidget::onAddPressed);

    layout->addLayout(header);

    this->m_tabwidget = new QTabWidget(this);
    this->m_tabwidget->setTabsClosable(true);

    connect(this->m_tabwidget, &QTabWidget::tabCloseRequested,
            this, &IntegerTabWidget::onTabClosed);

    this->m_key = 1;
    for (auto param = map->parameters.begin(); param != map->parameters.end(); ++param)
    {
        this->addParam(*param);

        ECFObject* obj = static_cast<ECFObject*>(*param);
        int k = QString::fromStdString(obj->name).toInt();
        if (k > this->m_key) this->m_key = k + 1;
    }

    layout->addWidget(this->m_tabwidget);
}

void IntegerTabWidget::onTabClosed(int index)
{
    QString key = this->m_tabwidget->tabText(index);
    this->m_map->dropParameter(
                    this->m_map->getParameter(key.toStdString())
                );
    this->m_tabwidget->removeTab(index);
}

void IntegerTabWidget::onAddPressed()
{
    QString key = QString::number(this->m_key++);

    this->addParam(this->m_map->getParameter(key.toStdString()));

    this->m_tabwidget->setCurrentIndex(this->m_tabwidget->count() - 1);
}

void IntegerTabWidget::addParam(ECFParameter *param)
{
    ECFObject* obj = static_cast<ECFObject*>(param);
    ECFObjectWidget* w = this->m_factory->create(obj);
    w->init();
    QString label = QString::fromStdString(obj->name);
    this->m_tabwidget->addTab(w, label);
}
