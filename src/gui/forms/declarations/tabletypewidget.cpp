#include "tabletypewidget.h"

#include "../validators/validatordelegate.h"
#include "../validators/validatorfactory.h"

TableTypeWidget::TableTypeWidget(QWidget* parent) :
    TableWidget(2, TableType::headlines(), parent)
{
    this->mTable->setItemDelegateForColumn(0, new ValidatorDelegate(new NumberValidatorFactory()));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new NumberValidatorFactory()));
}
