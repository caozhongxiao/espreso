#include "plot.h"
#include "ui_plot.h"
#include <QtMath>
#include <QGLFormat>
#include <QGraphicsTextItem>

Plot::Plot(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Plot)
{
    ui->setupUi(this);

    this->scene = new QGraphicsScene(this);
    scene->setSceneRect(0, 0, this->width() - 5, this->height() - 5);
    ui->view->setScene(scene);
    //ui->view->setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));

    fnXLeftBoundary = -20;
    fnXRightBoundary = 20;
    fnYTopBoundary = 1.5;
    fnYBottomBoundary = -1.5;
    fnXAxisLen = qAbs(fnXRightBoundary - fnXLeftBoundary);
    fnYAxisLen = qAbs(fnYTopBoundary - fnYBottomBoundary);
    sceneFnXRatio = this->width() / fnXAxisLen;
    sceneFnYRatio = this->height() / fnYAxisLen;
    qreal fnXAxisPointNum = 10;
    qreal fnYAxisPointNum = 10;

    this->font = QFont();
    this->fontSize = 0;
    if (font.pixelSize() != -1)
        this->fontSize = font.pixelSize();
    else if (font.pointSize() != -1)
        this->fontSize = font.pointSize();
    else if (font.pointSizeF() != -1)
        this->fontSize = font.pixelSize();
    else
        qWarning("%s", tr("Plot: Cannot retrieve the font size!").toStdString().c_str());

    this->fontSizeHalf = this->fontSize / 2;

    // Axises
    scene->addLine(this->width() / 2, 0, this->width() / 2, this->height());
    scene->addLine(0, this->height() / 2, this->width(), this->height() / 2);

    // Axis labels
    this->drawXAxisLabels(fnXAxisPointNum);
    this->drawYAxisLabels(fnYAxisPointNum);

    // Plot function
    qreal inc = qAbs(fnXLeftBoundary) / 100000;
    for (qreal x = fnXLeftBoundary; x < fnXRightBoundary; x += inc)
    {
        qreal y = this->fn(x);
        qreal sceneX = this->fnXToScene(x);
        qreal sceneY = this->fnYToScene(y);
        this->drawPoint(QPointF(sceneX, sceneY));
    }
}

Plot::~Plot()
{
    delete ui;
}

qreal Plot::fnXToScene(qreal x)
{
    return (x - this->fnXLeftBoundary) * this->sceneFnXRatio;
}

qreal Plot::fnYToScene(qreal y)
{
    return this->height() - (y - this->fnYBottomBoundary) * this->sceneFnYRatio;
}

qreal Plot::fn(qreal x)
{
    //return x;
    //return qPow(x, 2);
    return qSin(x);
    //return qSin(x) / x;
    //return qLn(x);
}

void Plot::drawPoint(QPointF p)
{
    if (p.rx() < 0 || p.rx() > this->width())
        return;
    if (p.ry() < 0 || p.ry() > this->height())
        return;

    this->scene->addRect(p.rx(), p.ry(), 0.5, 0.5, QPen(Qt::blue));
}

void Plot::drawXAxisLabels(int labelsCount, int labelPointLength)
{
    qreal labelPointRadius = labelPointLength / 2;
    qreal shift = this->fnXAxisLen / labelsCount;
    qreal y0 = this->fnYToScene(0);
    qreal y0Shifted = y0 + this->fontSizeHalf / 2;
    for (qreal x = this->fnXLeftBoundary + shift; x < this->fnXRightBoundary; x += shift)
    {
        qreal rounded = qRound(x * 1000) / 1000.0;
        if (rounded == 0) continue;
        QString content = QString::number(rounded);
        QGraphicsTextItem* text = scene->addText(content, this->font);
        qreal sceneX = this->fnXToScene(x);
        qreal textShift = content.length() * this->fontSizeHalf;
        text->setPos(sceneX - textShift, y0Shifted);
        scene->addLine(sceneX, y0 - labelPointRadius, sceneX, y0 + labelPointRadius);
    }
}

void Plot::drawYAxisLabels(int labelsCount, int labelPointLength)
{
    qreal labelPointRadius = labelPointLength / 2;
    qreal shift = this->fnYAxisLen / labelsCount;
    qreal x0 = this->fnXToScene(0);
    for (qreal y = this->fnYBottomBoundary + shift; y < this->fnYTopBoundary; y += shift)
    {
        qreal rounded = qRound(y * 1000) / 1000.0;
        if (rounded == 0) continue;
        QString content = QString::number(rounded);
        QGraphicsTextItem* text = scene->addText(content, this->font);
        qreal sceneY = this->fnYToScene(y);
        qreal len = content.length();
        if (content.at(0) == '-') len -= 0.5;
        text->setPos(x0 - len * this->fontSize, sceneY - this->fontSize);
        scene->addLine(x0 - labelPointRadius, sceneY, x0 + labelPointRadius, sceneY);
    }
}
