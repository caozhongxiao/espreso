#ifndef ECFOBJECTWIDGET_H
#define ECFOBJECTWIDGET_H

#include <QWidget>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "isavableobject.h"
#include "ivalidatableobject.h"
#include "optionhandler.h"
#include "boolhandler.h"
#include "../declarations/material/materialpropertytablewidget.h"
#include "formwidget.h"
#include "regionpropertywidget.h"

namespace espreso
{

namespace Ui {
class ECFObjectWidget;
}

class ECFObjectWidget : public QWidget,
        public IValidatableObject, public ISavableObject
{
    Q_OBJECT

public:
    explicit ECFObjectWidget(ECFObject* object, QWidget *parent = 0);
    virtual ~ECFObjectWidget();
    void init();

    virtual bool isValid() override;
    virtual QString errorMessage() override;

    virtual void save() override;

    void setDrawHeadline(bool draw);

protected slots:
    void redraw();

protected:
    Ui::ECFObjectWidget *ui;

    ECFObject* m_obj;

    QWidget* m_container;
    QWidget* m_widget;

    QVector<ISavableObject*> m_savables;
    QVector<IValidatableObject*> m_validatables;
    bool validate();

    bool m_draw_headlines = true;

    QString m_errormsg;

    virtual QWidget* initContainerWidget();
    virtual QWidget* initContainer() = 0;

    virtual void drawObject(ECFObject*);
    void drawMe();

    void processParameters(ECFObject*, QWidget*);
    virtual FormWidget* processOptionEnum(ECFParameter*, FormWidget*, QWidget*);
    virtual FormWidget* processBool(ECFParameter*, FormWidget*, QWidget*);
    virtual MaterialPropertyTableWidget* processExpression(ECFParameter*, MaterialPropertyTableWidget*, QWidget*);
    virtual FormWidget* processString(ECFParameter*, FormWidget*, QWidget*);
    virtual FormWidget* processFloat(ECFParameter*, FormWidget*, QWidget*);
    virtual FormWidget* processNonnegativeInteger(ECFParameter*, FormWidget*, QWidget*);
    virtual FormWidget* processPositiveInteger(ECFParameter*, FormWidget*, QWidget*);
    virtual FormWidget* processRegion(ECFParameter*, FormWidget*, QWidget*);

    OptionHandler* createOption(ECFParameter*, QWidget* = 0, bool = true);
    BoolHandler* createBool(ECFParameter*, QWidget* = 0);
    FormWidget* createFormWidget(QWidget*, FormWidget*);
    void createHeadline(ECFObject*, QWidget*);
};

}

#endif // ECFOBJECTWIDGET_H
