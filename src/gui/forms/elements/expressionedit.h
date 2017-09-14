#ifndef EXPRESSIONEDIT_H
#define EXPRESSIONEDIT_H

#include <QLineEdit>
#include <QFocusEvent>

class ExpressionEdit : public QLineEdit
{
    Q_OBJECT

public:
    ExpressionEdit(QWidget* parent = nullptr);
    ExpressionEdit(const QString& contents, QWidget* parent = nullptr);

    bool isValid();

    static bool validate(const QString& expr);

protected:
    virtual void focusOutEvent(QFocusEvent* e) override;

private:
    bool mValidState;
};

#endif // EXPRESSIONEDIT_H
