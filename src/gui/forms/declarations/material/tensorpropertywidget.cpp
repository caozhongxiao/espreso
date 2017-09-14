#include "tensorpropertywidget.h"
#include "ui_tensorpropertywidget.h"

using namespace espreso;

TensorPropertyWidget::TensorPropertyWidget(ECFObject* property, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TensorPropertyWidget)
{
    ui->setupUi(this);

    this->mProperty = property;

    connect(ui->cmbModel, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &TensorPropertyWidget::onIndexChanged);

    ui->lblName->setText(QString::fromStdString(property->metadata.description.at(0)));

    this->mWidget = new MaterialPropertyTableWidget(this);
    ui->layout->addWidget(mWidget);

    for (auto parameter = property->parameters.cbegin();
         parameter != property->parameters.cend()
         && (*parameter)->metadata.isallowed();
         ++parameter)
    {
        if ((*parameter)->metadata.datatype.at(0) == ECFDataType::OPTION)
        {
            this->mModel = (*parameter);
            for (auto option = (*parameter)->metadata.options.cbegin();
                 option != (*parameter)->metadata.options.cend()
                 && option->isallowed();
                 ++option)
            {
                ui->cmbModel->addItem(QString::fromStdString(option->name));
                this->mOptions.append(option->name);

                if ((*parameter)->getValue().compare(option->name) == 0)
                    this->onIndexChanged(0);
            }

            break;
        }
    }

}

TensorPropertyWidget::~TensorPropertyWidget()
{
    delete ui;
}

void TensorPropertyWidget::onIndexChanged(int index)
{
    if (index >= mOptions.size())
        return;

    MaterialPropertyTableWidget* w = this->mWidget;
    ui->layout->removeWidget(w);
    w->setParent(nullptr);
    delete w;

    this->mWidget = new MaterialPropertyTableWidget(this);
    ui->layout->addWidget(mWidget);

    this->mModel->setValue(mOptions.at(index));

    for (auto parameter = mProperty->parameters.cbegin();
         parameter != mProperty->parameters.cend()
         && (*parameter)->metadata.isallowed();
         ++parameter)
    {
        if ((*parameter)->metadata.datatype.at(0) == ECFDataType::EXPRESSION)
        {
            this->mWidget->addProperty( *parameter );
        }
    }
}
