#ifndef PIECEWISETYPEWIDGET_H
#define PIECEWISETYPEWIDGET_H

#include "tablewidget.h"
#include "../../data/datatype.h"

class PiecewiseTypeWidget : public TableWidget
{
    Q_OBJECT

signals:
    void validStateChanged(bool valid);

private slots:
    void changeValidState(bool valid);

public:
    PiecewiseTypeWidget(QWidget* parent = 0);
};

#endif // PIECEWISETYPEWIDGET_H
