#ifndef REGIONOBJECTWIDGET_H
#define REGIONOBJECTWIDGET_H

#include "ecfobjectwidget.h"

namespace espreso
{

class RegionObjectWidget : public ECFObjectWidget
{
    Q_OBJECT
public:
    explicit RegionObjectWidget(ECFObject* obj, QWidget *parent = nullptr);

protected:
    virtual QWidget* initContainer() override;
    virtual void drawObject(ECFObject*) override;
};

}

#endif // REGIONOBJECTWIDGET_H
