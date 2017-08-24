#ifndef PIECEWISETYPEWIDGET_H
#define PIECEWISETYPEWIDGET_H

#include "tablewidget.h"
#include "../../data/datatype.h"

class PiecewiseTypeWidget : public TableWidget
{
    Q_OBJECT

public:
    PiecewiseTypeWidget(QWidget* parent = 0);
    bool isValid();
};

#endif // PIECEWISETYPEWIDGET_H
