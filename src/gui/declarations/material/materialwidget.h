#ifndef MATERIALWIDGET_H
#define MATERIALWIDGET_H

#include "../../elements/ecfobjectwidget.h"
#include "../../elements/formwidget.h"

namespace espreso
{

class MaterialWidget : public ECFObjectWidget
{
    Q_OBJECT
public:
    MaterialWidget(MaterialConfiguration* material,
                   const QVector<std::string>& materialNames,
                   QWidget *parent = 0);

    bool isValid() override;

protected:
    QWidget* initContainer() override;
    void drawObject(ECFObject*) override;

private:
    QVector<std::string> m_names;
    QWidget* m_widget;

    void processParameters(ECFObject*, QWidget*);
    void drawHeadline(ECFObject*, QWidget*);

    FormWidget* createFormWidget(QWidget*, FormWidget* = nullptr);
};

}
#endif // MATERIALWIDGET_H
