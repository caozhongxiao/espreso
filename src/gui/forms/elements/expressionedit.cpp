#include "expressionedit.h"

#include "../../../basis/expression/expression.h"
#include "../../data/common.h"

ExpressionEdit::ExpressionEdit(QWidget* parent) : QLineEdit(parent)
{
    this->setText("0");
    this->mValidState = true;
}

ExpressionEdit::ExpressionEdit(const QString& contents, QWidget* parent) :
    QLineEdit(contents, parent)
{
    this->setText("0");
    this->mValidState = true;
}

bool ExpressionEdit::validate(const QString& expr)
{
    return espreso::Expression::isValid(expr.toStdString(), Common::fnVariables());
}

bool ExpressionEdit::isValid()
{
    return this->mValidState;
}

void ExpressionEdit::focusOutEvent(QFocusEvent* e)
{
    bool valid = ExpressionEdit::validate(this->text());

    if (valid)
    {
        this->setStyleSheet("ExpressionEdit{color: black;}");
    }
    else
    {
        this->setStyleSheet("ExpressionEdit{color: red; font-weight: bold;}");
    }

    this->mValidState = valid;

    QLineEdit::focusOutEvent(e);
}
