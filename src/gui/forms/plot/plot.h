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

    qreal fnXLeftBoundary;
    qreal fnXRightBoundary;
    qreal fnYTopBoundary;
    qreal fnYBottomBoundary;
    qreal fnXAxisLen;
    qreal fnYAxisLen;
    qreal sceneFnXRatio;
    qreal sceneFnYRatio;

    qreal fnXToScene(qreal x);
    qreal fnYToScene(qreal y);

    qreal fn(qreal x);

    void drawPoint(QPointF p);
    void drawXAxisLabels(int labelsCount, int labelPointLength = 10, const QFont& font = QFont());
    void drawYAxisLabels(int labelsCount, int labelPointLength = 10, const QFont& font = QFont());
};

#endif // PLOT_H
