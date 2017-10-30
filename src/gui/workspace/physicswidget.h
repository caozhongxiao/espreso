#ifndef PHYSICSWIDGET_H
#define PHYSICSWIDGET_H

#include <QPushButton>
#include <QComboBox>

#include "../elements/ecfobjectwidget.h"
#include "../elements/fieldhandler.h"

namespace espreso
{

class PhysicsWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    explicit PhysicsWidget(ECFConfiguration* ecf, QWidget* parent = 0);

    ECFObject* activePhysics();

signals:
    void loadstepsChanged(int loadsteps);
    void physicsChanged(ECFObject* physics);

protected:
    QWidget* initContainer() override;
    void drawObject(ECFObject*) override;

private slots:
    void onPhysicsChange(int index);
    void onLoadstepsChange(int loadsteps);

private:
    ECFConfiguration* m_ecf;
    QWidget* m_widget;

    QComboBox* m_physics;

    void processParameters(ECFObject* obj, QWidget* widget);
    ECFObject* physics(int index);
};

}

#endif // PHYSICSWIDGET_H
