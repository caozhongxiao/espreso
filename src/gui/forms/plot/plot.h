#ifndef PLOT_H
#define PLOT_H

#include <QWidget>
#include <QGraphicsScene>

namespace Ui {
class Plot;
}

class Plot : public QWidget
{
    Q_OBJECT

public:
    explicit Plot(QWidget *parent = 0);
    ~Plot();

private:
    Ui::Plot *ui;

    QGraphicsScene* scene;
    qreal ratio;

    qreal coordinateToX(qreal coordinate);
    qreal yToCoordinate(qreal y);
    qreal fn(qreal x);
    void drawPoint(QPointF p);
};

#endif // PLOT_H
