#ifndef PLOT_H
#define PLOT_H

#include <QWidget>
#include <QGraphicsScene>

namespace Ui {
class Plot;
}

struct PlotLabel
{
    QString text;
    qreal width;
    qreal x;
    qreal y;
};

class Plot : public QWidget
{
    Q_OBJECT

public:
    explicit Plot(QWidget *parent = 0);
    ~Plot();

private:
    Ui::Plot *ui;

    QGraphicsScene* scene;

    QFont font;
    qreal fontSize;
    qreal fontSizeHalf;

    qreal fnXLeftBoundary;
    qreal fnXRightBoundary;
    qreal fnYTopBoundary;
    qreal fnYBottomBoundary;
    qreal fnXAxisLen;
    qreal fnYAxisLen;
    qreal rectFnXRatio;
    qreal rectFnYRatio;

    QGraphicsRectItem* rect;
    qreal rectX;
    qreal rectY;
    qreal rectWidth;
    qreal rectHeight;
    qreal plotXMiddle;
    qreal plotYMiddle;
    qreal mainRectX();

    qreal xAxisPrecision;
    qreal yAxisPrecision;
    qreal computePrecision(qreal intervalLength);

    QVector<PlotLabel> yLabels;
    QVector<PlotLabel> xLabels;

    qreal fnXToRect(qreal x);
    qreal fnYToRect(qreal y);
    qreal rectXToFn(qreal x);
    qreal rectYToFn(qreal y);

    qreal fn(qreal x);

    void drawPoint(QPointF p);
    void drawAxises();
    void drawFn();
    void drawXLabels(int labelPointLength = 10);
    void drawYLabels(int labelPointLength = 10);
};

#endif // PLOT_H
