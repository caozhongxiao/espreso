#ifndef MATRIXPROPERTYWIDGET_H
#define MATRIXPROPERTYWIDGET_H

#include <QWidget>

namespace Ui {
class MatrixPropertyWidget;
}

class MatrixPropertyWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MatrixPropertyWidget(QWidget *parent = 0);
    ~MatrixPropertyWidget();

private:
    Ui::MatrixPropertyWidget *ui;
};

#endif // MATRIXPROPERTYWIDGET_H
