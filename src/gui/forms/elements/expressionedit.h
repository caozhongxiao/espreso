#ifndef EXPRESSIONEDIT_H
#define EXPRESSIONEDIT_H

#include <QLineEdit>
#include <QFocusEvent>
#include "ivalidatableobject.h"

namespace espreso
{

    class ExpressionEdit : public QLineEdit, public IValidatableObject
    {
        Q_OBJECT

    public:
        ExpressionEdit(QWidget* parent = nullptr);
        ExpressionEdit(const QString& contents, QWidget* parent = nullptr);

        bool isValid() override;

        static bool validate(const QString& expr);

    protected:
        virtual void focusOutEvent(QFocusEvent* e) override;

    private:
        bool mValidState;
    };

}

#endif // EXPRESSIONEDIT_H
