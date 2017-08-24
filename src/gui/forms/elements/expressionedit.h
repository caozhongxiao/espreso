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

signals:
    void validStateChanged(bool isValid);

protected:
    virtual void focusOutEvent(QFocusEvent* e) override;

private:
    bool mValidState;
};

#endif // EXPRESSIONEDIT_H
