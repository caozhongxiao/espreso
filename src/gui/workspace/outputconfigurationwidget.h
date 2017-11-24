#ifndef OUTPUTCONFIGURATIONWIDGET_H
#define OUTPUTCONFIGURATIONWIDGET_H

#include "../elements/ecfobjectwidget.h"

namespace espreso
{

class OutputConfigurationWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    OutputConfigurationWidget(ECFObject* output, QWidget* parent = 0);

protected:
    virtual QWidget* initContainer() override;
    virtual void drawObject(ECFObject*) override;

private:
    QWidget* m_widget;
};

}

#endif // OUTPUTCONFIGURATIONWIDGET_H
