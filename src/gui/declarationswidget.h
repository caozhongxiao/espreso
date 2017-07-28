#ifndef DECLARATIONSWIDGET_H
#define DECLARATIONSWIDGET_H

#include <QWidget>

namespace Ui {
class DeclarationsWidget;
}

class DeclarationsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DeclarationsWidget(QWidget *parent = 0);
    ~DeclarationsWidget();

private:
    Ui::DeclarationsWidget *ui;
};

#endif // DECLARATIONSWIDGET_H
