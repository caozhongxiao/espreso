#include "piecewisetypewidget.h"

#include "../elements/expressioneditdelegate.h"
#include "../validators/validatordelegate.h"

PiecewiseTypeWidget::PiecewiseTypeWidget(QWidget* parent) :
    TableWidget(3, PiecewiseFunctionType::headlines(), parent)
{
    this->mTable->setItemDelegateForColumn(0, new ValidatorDelegate(new NumberValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new NumberValidatorFactory(), this));

    ExpressionEditDelegate* eed = new ExpressionEditDelegate(this);
    connect(eed, &ExpressionEditDelegate::validStateChanged, this, &PiecewiseTypeWidget::changeValidState);
    this->mTable->setItemDelegateForColumn(2, eed);
}

void PiecewiseTypeWidget::changeValidState(bool valid)
{
    emit validStateChanged(valid);
}
