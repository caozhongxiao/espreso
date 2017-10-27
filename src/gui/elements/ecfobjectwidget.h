#ifndef ECFOBJECTWIDGET_H
#define ECFOBJECTWIDGET_H

#include <QWidget>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "isavableobject.h"
#include "ivalidatableobject.h"
#include "optionhandler.h"

namespace espreso
{

namespace Ui {
class ECFObjectWidget;
}

class ECFObjectWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ECFObjectWidget(ECFObject* object, QWidget *parent = 0);
    virtual ~ECFObjectWidget();
    void init();

protected slots:
    void redraw();

protected:
    Ui::ECFObjectWidget *ui;

    ECFObject* m_obj;

    QWidget* m_container;

    QVector<ISavableObject*> m_savables;
    QVector<IValidatableObject*> m_validatables;
    bool validate();

    virtual QWidget* initContainer() = 0;

    virtual void drawObject(ECFObject*) = 0;
    void drawMe();

    OptionHandler* createOption(ECFParameter*, QWidget* = 0, bool = true);
};

}

#endif // ECFOBJECTWIDGET_H
