#ifndef PHYSICSWIDGET_H
#define PHYSICSWIDGET_H

#include "../elements/ecfobjectwidget.h"

namespace espreso
{

class PhysicsWidget : public ECFObjectWidget
{
public:
    PhysicsWidget(ECFConfiguration* ecf, QWidget* parent = 0);

protected:
    QWidget* initContainer() override;
    void drawObject(ECFObject*) override;

private slots:
    void onPhysicsChange(int index);

private:
    ECFConfiguration* m_ecf;
    QWidget* m_widget;

    void processParameters(ECFObject* obj, QWidget* widget);
    ECFObject* physics(int index);
};

}

#endif // PHYSICSWIDGET_H
