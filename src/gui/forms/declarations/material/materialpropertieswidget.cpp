#include "materialpropertieswidget.h"
#include "ui_materialpropertieswidget.h"

MaterialPropertiesWidget::MaterialPropertiesWidget(const QVector<TensorProperty>& tensors,
                                                   const QVector<ScalarProperty>& scalars,
                                                   QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MaterialPropertiesWidget)
{
    ui->setupUi(this);

    // TENSORS
    foreach (TensorProperty tp, tensors) {
//        TensorPropertyWidget* w = new TensorPropertyWidget(tp, this);
//        this->tensorWidgets.append(w);
//        ui->layoutTensor->addWidget(w);
    }

    // SCALARS
    this->scalarWidget = new MaterialPropertyTableWidget(this);
    foreach (ScalarProperty sp, scalars) {
        //scalarWidget->addProperty(&sp);
    }
    ui->layoutScalar->addWidget(scalarWidget);
}

MaterialPropertiesWidget::~MaterialPropertiesWidget()
{
    delete ui;
}
