#ifndef FIXEDECFOBJECTWIDGET_H
#define FIXEDECFOBJECTWIDGET_H

#include "ecfobjectwidget.h"

namespace espreso
{

class FixedECFObjectWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    FixedECFObjectWidget(ECFObject* obj, QWidget* parent = 0);

    void setDrawHeadline(bool draw);

protected:
    virtual QWidget* initContainer() override;
    virtual void drawObject(ECFObject*) override;

private:
    bool m_draw_headlines;
};

}

#endif // FIXEDECFOBJECTWIDGET_H
