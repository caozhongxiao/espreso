#ifndef LOADSTEPWIDGET_H
#define LOADSTEPWIDGET_H

#include "../elements/ecfobjectwidget.h"
#include "../elements/fieldhandler.h"

namespace espreso
{

class LoadstepWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    explicit LoadstepWidget(size_t id, ECFObject* physics, QWidget* parent = 0);

protected:
    QWidget* initContainer() override;
    void drawObject(ECFObject*) override;

private:
    ECFObject* m_physics;
    ECFParameter* m_loadstep;

    QWidget* m_widget;

    void processParameters(ECFObject* obj, QWidget* widget);
};

}

#endif // LOADSTEPWIDGET_H
