#include "plot.h"
#include "ui_plot.h"
#include <QtMath>

Plot::Plot(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Plot)
{
    ui->setupUi(this);

    this->scene = new QGraphicsScene(this);
    scene->setSceneRect(0, 0, this->width() - 5, this->height() - 5);
    ui->view->setScene(scene);

    qreal fnXLeftBoundary = -5;
    qreal fnXRightBoundary = 5;
    qreal fnYTopBoundary = 5;
    qreal fnYBottomBoundary = -5;
    qreal fnXAxisLen = qAbs(fnXRightBoundary - fnXLeftBoundary);
    qreal fnYAxisLen = qAbs(fnYTopBoundary - fnYBottomBoundary);
    qreal sceneFnXRatio = this->width() / fnXAxisLen;
    qreal sceneFnYRatio = this->height() / fnYAxisLen;


    this->ratio = 10;

    // Axises
    scene->addLine(this->width() / 2, 0, this->width() / 2, this->height());
    scene->addLine(0, this->height() / 2, this->width(), this->height() / 2);

    qreal inc = qAbs(fnXLeftBoundary) / 100;
    for (qreal x = fnXLeftBoundary; x < fnXRightBoundary; x += inc)
    {
        qreal y = this->fn(x);
        qreal sceneX = (x - fnXLeftBoundary) * sceneFnXRatio;
        qreal sceneY = this->height() - (y - fnYBottomBoundary) * sceneFnYRatio;
        qInfo("x: %f, y: %f / sX: %f, sY: %f", x, y, sceneX, sceneY);
        this->drawPoint(QPointF(sceneX, sceneY));
    }
}

Plot::~Plot()
{
    delete ui;
}

qreal Plot::coordinateToX(qreal coordinate)
{
    return coordinate - (this->width() / 2);
}

qreal Plot::yToCoordinate(qreal y)
{
    return this->height() / 2 - y;
}

qreal Plot::fn(qreal x)
{
    //return x;
    return qPow(x, 2);
    //return qCos(x);
    //return qLn(x);
}

void Plot::drawPoint(QPointF p)
{
    if (p.rx() < 0 || p.rx() > this->width())
        return;
    if (p.ry() < 0 || p.ry() > this->height())
        return;

    this->scene->addRect(p.rx(), p.ry(), 0.5, 0.5);
}
