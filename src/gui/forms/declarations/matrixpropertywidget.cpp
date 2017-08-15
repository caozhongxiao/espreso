#include "matrixpropertywidget.h"
#include "ui_matrixpropertywidget.h"

MatrixPropertyWidget::MatrixPropertyWidget(const QHash<QString, Variable>& varDict, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MatrixPropertyWidget)
{
    ui->setupUi(this);

    this->varDict = varDict;

    ui->frDetail->hide();
    ui->tblMatrix->setEditTriggers(QAbstractItemView::NoEditTriggers);

    this->matrixModel = new QStandardItemModel(3, 3);
    ui->tblMatrix->setModel(this->matrixModel);
    QStringList tblHeader;
    tblHeader << "x" << "y" << "z";
    this->matrixModel->setHorizontalHeaderLabels(tblHeader);
    QStringList tblSide;
    tblSide << "x" << "y" << "z";
    this->matrixModel->setVerticalHeaderLabels(tblSide);

    this->matrixItemWidget = new DataTypeWidget(this->varDict, ui->frEdit);
    ui->frEdit->layout()->addWidget(this->matrixItemWidget);

    QVector<DataType*> mRow;
    mRow.append(new ConstantType("0"));
    mRow.append(new ConstantType("0"));
    mRow.append(new ConstantType("0"));
    this->matrix.append(mRow);
    this->matrix.append(mRow);
    this->matrix.append(mRow);
}

MatrixPropertyWidget::MatrixPropertyWidget(MaterialProperty* property, const QHash<QString, Variable>& varDict, QWidget* parent) :
    MatrixPropertyWidget(varDict, parent)
{
    this->setData(property);
}

MatrixPropertyWidget::~MatrixPropertyWidget()
{
    delete ui;
}

void MatrixPropertyWidget::setData(MaterialProperty* property)
{
    ui->frDetail->hide();
    this->mProperty = property;
    this->matrixModel->setRowCount(0);
    this->matrixModel->setRowCount(3);
    property->accept(this);
}

void MatrixPropertyWidget::visit(BasicProperty& p) {}

void MatrixPropertyWidget::visit(IsotropicProperty& p)
{
    this->propertyType = 1;
    QPushButton* btn = new QPushButton(tr("Edit"), ui->tblMatrix);
    connect(btn, &QPushButton::pressed,
            [=] () { matrixBtnPressed(0, 0); });
    QModelIndex btnIndex = this->matrixModel->index(0, 0);
    QModelIndex yy = this->matrixModel->index(1,1);
    QModelIndex zz = this->matrixModel->index(2,2);
    this->matrixModel->setData(yy, QVariant(tr("Inherited from xx")));
    this->matrixModel->setData(zz, QVariant(tr("Inherited from xx")));
    ui->tblMatrix->setIndexWidget(btnIndex, btn);
    this->matrix[0][0] = p.model().kxx;
}

void MatrixPropertyWidget::visit(DiagonalProperty& p)
{
    this->propertyType = 2;
    QPushButton* btn1 = new QPushButton(tr("Edit"), ui->tblMatrix);
    QPushButton* btn2 = new QPushButton(tr("Edit"), ui->tblMatrix);
    QPushButton* btn3 = new QPushButton(tr("Edit"), ui->tblMatrix);
    connect(btn1, &QPushButton::pressed,
            [=] () { matrixBtnPressed(0, 0); });
    connect(btn2, &QPushButton::pressed,
            [=] () { matrixBtnPressed(1, 1); });
    connect(btn3, &QPushButton::pressed,
            [=] () { matrixBtnPressed(2, 2); });
    QModelIndex xx = this->matrixModel->index(0, 0);
    QModelIndex yy = this->matrixModel->index(1, 1);
    QModelIndex zz = this->matrixModel->index(2, 2);
    ui->tblMatrix->setIndexWidget(xx, btn1);
    ui->tblMatrix->setIndexWidget(yy, btn2);
    ui->tblMatrix->setIndexWidget(zz, btn3);
    this->matrix[0][0] = p.model().kxx;
    this->matrix[1][1] = p.model().kyy;
    this->matrix[2][2] = p.model().kzz;
}

void MatrixPropertyWidget::visit(SymmetricProperty& p)
{
    this->propertyType = 3;
    QVector<QVector<int> > items(3);
    items[0].append(0);
    items[0].append(1);
    items[0].append(2);
    items[1].append(1);
    items[1].append(2);
    items[2].append(2);
    for (int i = 0; i < 3; i++)
    {
        QVector<int> row = items.at(i);
        for (auto j = row.cbegin(); j != row.cend(); ++j)
        {
            QPushButton* btn = new QPushButton(tr("Edit"), ui->tblMatrix);
            int col = *j;
            connect(btn, &QPushButton::pressed,
                    [=] () { matrixBtnPressed(i, col); });
            QModelIndex index = this->matrixModel->index(i, *j);
            ui->tblMatrix->setIndexWidget(index, btn);
            this->matrix[i][*j] = p.modelData().at( i * 3 + (*j) );
        }
    }
}

void MatrixPropertyWidget::visit(AnisotropicProperty& p)
{
    this->propertyType = 4;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            QPushButton* btn = new QPushButton(tr("Edit"), ui->tblMatrix);
            connect(btn, &QPushButton::pressed,
                    [=] () { matrixBtnPressed(i, j); });
            QModelIndex index = this->matrixModel->index(i, j);
            ui->tblMatrix->setIndexWidget(index, btn);
            this->matrix[i][j] = p.modelData().at( i * 3 + j );
        }
    }
}

void MatrixPropertyWidget::matrixBtnPressed(int row, int column)
{
    this->pRow = row;
    this->pColumn = column;
    this->matrixItemWidget->setDataType(this->matrix[row][column]);
    ui->frDetail->show();
}

void MatrixPropertyWidget::on_btnSave_pressed()
{
    DataType* result;

    if ((result = this->matrixItemWidget->checkedData()) == nullptr)
        return;

    this->matrix[pRow][pColumn] = result;
    ui->frDetail->hide();
}

MaterialProperty* MatrixPropertyWidget::data()
{
    if (this->propertyType == 1)
    {
        IsotropicProperty* ip = (IsotropicProperty*)this->mProperty;
        ip->model().kxx = this->matrix[0][0];
    }
    else if (this->propertyType == 2)
    {
        DiagonalProperty* dp = (DiagonalProperty*)this->mProperty;
        dp->model().kxx = this->matrix[0][0];
        dp->model().kyy = this->matrix[1][1];
        dp->model().kzz = this->matrix[2][2];
    }
    else if (this->propertyType == 3)
    {
        SymmetricProperty* sp = (SymmetricProperty*)this->mProperty;
        sp->model().kxx = this->matrix[0][0];
        sp->model().kyy = this->matrix[1][1];
        sp->model().kzz = this->matrix[2][2];
        sp->model().kxy = this->matrix[0][1];
        sp->model().kxz = this->matrix[0][2];
        sp->model().kyz = this->matrix[2][2];
    }
    else if (this->propertyType == 4)
    {
        AnisotropicProperty* ap = (AnisotropicProperty*)this->mProperty;
        ap->model().kxx = this->matrix[0][0];
        ap->model().kyy = this->matrix[1][1];
        ap->model().kzz = this->matrix[2][2];
        ap->model().kxy = this->matrix[0][1];
        ap->model().kxz = this->matrix[0][2];
        ap->model().kyz = this->matrix[2][2];
        ap->model().kyx = this->matrix[1][0];
        ap->model().kzx = this->matrix[2][0];
        ap->model().kzy = this->matrix[2][1];
    }
    else
    {
        qWarning("%s", tr("Unknown or not supported type of MaterialProperty!").toStdString().c_str());
    }

    return this->mProperty;
}
