#ifndef INPUTWIDGET_H
#define INPUTWIDGET_H

#include "../elements/fixedecfobjectwidget.h"

namespace espreso
{

class InputWidget : public FixedECFObjectWidget
{
public:
    InputWidget(ECFObject* obj, QWidget* parent = 0);

protected:
    void drawObject(ECFObject*) override;
};

}

#endif // INPUTWIDGET_H
