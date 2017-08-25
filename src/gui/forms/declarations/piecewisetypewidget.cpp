#include "piecewisetypewidget.h"

#include "../elements/expressioneditdelegate.h"
#include "../elements/expressionedit.h"
#include "../validators/validatordelegate.h"

PiecewiseTypeWidget::PiecewiseTypeWidget(QWidget* parent) :
    TableWidget(3, PiecewiseFunctionType::headlines(), parent)
{
    this->mTable->setItemDelegateForColumn(0, new ValidatorDelegate(new NumberValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new NumberValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(2, new ExpressionEditDelegate(this));
}

bool PiecewiseTypeWidget::isValid()
{
    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QModelIndex index = this->mModel->index(row, 2);
        QString expression = index.data().toString();
        if (!ExpressionEdit::validate(expression))
            return false;
    }

    return TableWidget::isValid();
}
