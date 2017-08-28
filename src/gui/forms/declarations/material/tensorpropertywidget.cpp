#include "tensorpropertywidget.h"
#include "ui_tensorpropertywidget.h"

TensorPropertyWidget::TensorPropertyWidget(const TensorProperty& property, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TensorPropertyWidget)
{
    ui->setupUi(this);

    this->mProperty = property;
    ui->lblName->setText(mProperty.name());
    int index = 0;
    for (auto model = mProperty.modelBegin(); model != mProperty.modelEnd(); ++model)
    {
        ui->cmbModel->addItem(model->name());

        MaterialPropertyTableWidget* w = new MaterialPropertyTableWidget(this);

        foreach (TensorPropertyModelItem item, model->items()) {
            w->addProperty(&item);
        }

        if (mProperty.activeModel() != index)
            w->hide();

        ui->layout->addWidget(w);

        this->mWidgets.append(w);

        index++;
    }

    ui->cmbModel->setCurrentIndex(mProperty.activeModel());
    connect(ui->cmbModel, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &TensorPropertyWidget::onIndexChanged);
}

TensorPropertyWidget::~TensorPropertyWidget()
{
    this->mWidgets.clear();
    delete ui;
}

void TensorPropertyWidget::onIndexChanged(int index)
{
    foreach (MaterialPropertyTableWidget* w, mWidgets) {
        w->hide();
    }
    this->mWidgets.at(index)->show();
}
