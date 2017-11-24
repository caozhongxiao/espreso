#ifndef INTEGERTABWIDGET_H
#define INTEGERTABWIDGET_H

#include <QWidget>
#include <QTabWidget>

#include <memory>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "ecfobjectwidgetfactory.h"

namespace espreso
{

class IntegerTabWidget : public QWidget
{
    Q_OBJECT

public:
    IntegerTabWidget(ECFObject* map, std::unique_ptr<ECFObjectWidgetFactory> factory, QWidget* parent = 0);

private slots:
    void onTabClosed(int index);
    void onAddPressed();

private:
    QTabWidget* m_tabwidget;

    ECFObject* m_map;
    std::unique_ptr<ECFObjectWidgetFactory> m_factory;

    int m_key;

    void addParam(ECFParameter*);
};

}

#endif // INTEGERTABWIDGET_H
