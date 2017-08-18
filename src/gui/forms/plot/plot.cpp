#include "plot.h"
#include "ui_plot.h"
#include <QtMath>
#include <QGLFormat>
#include <QGraphicsTextItem>
#include <cmath>

Plot::Plot(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Plot)
{
    ui->setupUi(this);

    this->scene = new QGraphicsScene(this);
    scene->setSceneRect(0, 0, this->width() - 5, this->height() - 5);
    ui->view->setScene(scene);
    //ui->view->setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));

    this->fnXLeftBoundary = -5;
    this->fnXRightBoundary = 7;
    this->fnYTopBoundary = 5;
    this->fnYBottomBoundary = -3;
    this->fnXAxisLen = qAbs(this->fnXRightBoundary - this->fnXLeftBoundary);
    this->fnYAxisLen = qAbs(this->fnYTopBoundary - this->fnYBottomBoundary);

    this->xAxisPrecision = this->computePrecision(fnXAxisLen);
    this->yAxisPrecision = this->computePrecision(fnYAxisLen);

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

    this->rectY = this->fontSize * 3;
    this->rectHeight = this->height() - this->rectY - this->fontSize * 3;
    this->rectX = this->mainRectX() + this->fontSizeHalf;
    this->rectWidth = this->width() - 2 * this->rectX;
    this->rect = this->scene->addRect(rectX, rectY, rectWidth, rectHeight);
    this->rectFnXRatio = this->rectWidth / this->fnXAxisLen;
    this->rectFnYRatio = this->rectHeight / this->fnYAxisLen;

    // Labels
    this->drawXLabels();
    this->drawYLabels();

    this->drawAxises();
    this->drawFn();
}

Plot::~Plot()
{
    delete ui;
}


void Plot::drawAxises()
{
    // Axises
    QPen axisColor = QPen(Qt::gray);
    if ( (0 >= fnXLeftBoundary) && (0 <= fnXRightBoundary) )
    {
        qreal x = this->fnXToRect(0);
        scene->addLine(x, this->rectY,
                       x, this->rectHeight + this->rectY,
                       axisColor);
    }
    if ( (0 >= fnYBottomBoundary) && (0 <= fnYTopBoundary) )
    {
        qreal y = this->fnYToRect(0);
        scene->addLine(this->rectX, y,
                       this->rectWidth + this->rectX, y,
                       axisColor);
    }
}

void Plot::drawFn()
{
    // Plot function
    qreal inc = qAbs(this->fnXLeftBoundary) / 100000;
    for (qreal x = this->fnXLeftBoundary; x < this->fnXRightBoundary; x += inc)
    {
        qreal y = this->fn(x);
        qreal sceneX = this->fnXToRect(x);
        qreal sceneY = this->fnYToRect(y);
        this->drawPoint(QPointF(sceneX, sceneY));
    }
}

qreal Plot::fnXToRect(qreal x)
{
    return (x - this->fnXLeftBoundary) * this->rectFnXRatio + this->rectX;
}

qreal Plot::fnYToRect(qreal y)
{
    return this->rectHeight - (y - this->fnYBottomBoundary) * this->rectFnYRatio + this->rectY;
}

qreal Plot::rectXToFn(qreal x)
{
    return (x - this->rectX) / this->rectFnXRatio + this->fnXLeftBoundary;
}

qreal Plot::rectYToFn(qreal y)
{
    return ( this->rectY + this->rectHeight - y ) / this->rectFnYRatio
            + this->fnYBottomBoundary;
}

qreal Plot::fn(qreal x)
{
    //return x;
    //return qPow(x, 2);
    //return qSin(x);
    //return qSin(x) / x;
    //return qLn(x);
    return qSin(x) / x + qLn(x);
}

void Plot::drawPoint(QPointF p)
{
    if ( p.rx() < this->rectX || p.rx() > (this->rectWidth + this->rectX) )
        return;
    if ( p.ry() < this->rectY || p.ry() > (this->rectHeight + this->rectY) )
        return;

    this->scene->addRect(p.rx(), p.ry(), 0.5, 0.5, QPen(Qt::blue));
}

qreal Plot::mainRectX()
{
    qreal labelCount = floor(this->rectHeight / (this->fontSize * 2));
    qreal shift = this->fnYAxisLen / labelCount;
    qreal max = 0;
    for (qreal y = this->fnYBottomBoundary; y < this->fnYTopBoundary; y += shift)
    {
        qreal rounded = floor(y * yAxisPrecision) / yAxisPrecision;

        QString content = QString::number(rounded);

        qreal len = content.length();

        if (content.at(0) == '-' && len > 2) len -= 0.5;
        if (len <= 2) len = 2.5;

        len *= this->fontSize;

        if (len > max) max = len;

        PlotLabel lbl;
        lbl.text = content;
        lbl.width = len;
        lbl.y = rounded;
        lbl.x = this->rectX;
        this->yLabels.append(lbl);
    }
    return max;
}

void Plot::drawXLabels(int labelPointLength)
{
    qreal labelPointRadius = labelPointLength / 2;

    qreal labelCount = log10(xAxisPrecision) +
            ceil( fmax( log10(qAbs(fnXLeftBoundary)),
                         log10(qAbs(fnXRightBoundary)) ) );

    qreal shift = this->fnXAxisLen / (labelCount * this->fontSize);

    qreal y0 = this->rectY + this->rectHeight;
    for (qreal x = this->fnXLeftBoundary; x <= this->fnXRightBoundary; x += shift)
    {
        qreal rounded = floor(x * xAxisPrecision) / xAxisPrecision;

        QString content = QString::number(rounded);
        QGraphicsTextItem* text = scene->addText(content, this->font);

        qreal sceneX = this->fnXToRect(rounded);
        qreal textShift = content.length() * this->fontSizeHalf;
        text->setPos(sceneX - textShift, y0 + this->fontSizeHalf);

        scene->addLine(sceneX, y0 - labelPointRadius, sceneX, y0 + labelPointRadius);
    }
}

void Plot::drawYLabels(int labelPointLength)
{
    qreal labelPointRadius = labelPointLength / 2;
    qreal x0 = this->rectX;
    foreach (PlotLabel pl, this->yLabels)
    {
        QGraphicsTextItem* text = this->scene->addText(pl.text, this->font);
        qreal sceneY = this->fnYToRect(pl.y);
        text->setPos(x0 - pl.width, sceneY - this->fontSize);
        scene->addLine(x0 - labelPointRadius, sceneY, x0 + labelPointRadius, sceneY);
    }
}

qreal Plot::computePrecision(qreal intervalLength)
{
    if (intervalLength >= 10)
        return 1;

    if (intervalLength < 10 && intervalLength >= 1)
        return 10;

    qreal ret = 10;
    QString number = QString::number(intervalLength);

    QChar prev('0');
    foreach (QChar ch, number) {
        if (prev != '0' || prev != '.')
            break;
        if (ch == '0')
            ret *= 10;
        prev = ch;
    }

    return ret;
}
