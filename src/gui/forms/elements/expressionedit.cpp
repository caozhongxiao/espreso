#include "expressionedit.h"

#include "../../expression.h"
#include "../../data/common.h"

ExpressionEdit::ExpressionEdit(QWidget* parent) : QLineEdit(parent)
{
    this->setText("0");
    this->mValidState = true;
}

ExpressionEdit::ExpressionEdit(const QString& contents, QWidget* parent) :
    QLineEdit(contents, parent)
{
    this->mValidState = true;
}

void ExpressionEdit::focusOutEvent(QFocusEvent* e)
{
    bool valid = Expression::isValid(this->text().toStdString(), Common::fnVariables());

    if (valid)
    {
        this->setStyleSheet("color: black;");
    }
    else
    {
        this->setStyleSheet("color: red; font-weight: bold;");
    }

    if ((mValidState == true && valid == false)
            || (mValidState == false && valid == true))
    {
        emit validStateChanged(valid);
        this->mValidState = valid;
    }

    QLineEdit::focusOutEvent(e);
}
