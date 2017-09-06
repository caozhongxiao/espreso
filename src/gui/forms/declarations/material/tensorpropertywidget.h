#ifndef TENSORPROPERTYWIDGET_H
#define TENSORPROPERTYWIDGET_H

#include <QWidget>
#include "materialpropertytablewidget.h"
#include "../../../../config/configuration.h"

namespace Ui {
class TensorPropertyWidget;
}

class TensorPropertyWidget : public QWidget
{
    Q_OBJECT

public:
    explicit TensorPropertyWidget(ECFObject* property, QWidget *parent = 0);
    ~TensorPropertyWidget();

private slots:
    void onIndexChanged(int index);

private:
    Ui::TensorPropertyWidget *ui;

    MaterialPropertyTableWidget* mWidget;
    QVector<std::string> mOptions;
    ECFObject* mProperty;
    ECFParameter* mModel;
};

#endif // TENSORPROPERTYWIDGET_H
