#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include <QDialog>
#include <QStandardItemModel>
#include "../data/material.h"
#include "../data/variable.h"
#include "datatypewidget.h"
#include "matrixpropertywidget.h"

namespace Ui {
class MaterialDialog;
}

class MaterialDialog : public QDialog, public MaterialPropertyVisitor
{
    Q_OBJECT

public:
    explicit MaterialDialog(const QHash<QString, Material>& matDict, const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    MaterialDialog(const Material& material, const QHash<QString, Material>& matDict, const QHash<QString, Variable>& varDict, QWidget *parent = 0);
    ~MaterialDialog();

    Material data() const;

    void visit(BasicProperty& p) override;
    void visit(IsotropicProperty& p) override;
    void visit(DiagonalProperty& p) override;
    void visit(SymmetricProperty& p) override;
    void visit(AnisotropicProperty& p) override;

private slots:
    void on_btnPropAdd_pressed();

    void on_btnPropDel_pressed();

    void on_btnPropSave_pressed();

private:
    Ui::MaterialDialog *ui;

    QHash<QString, Material> matDict;
    QHash<QString, Variable> varDict;

    QVector<MaterialProperty*> properties;
    QHash<int, int> deletedProperties;

    int activeProperty = 0;
    int activePropertyWidget = 0;
    int serveBasicPropertyWidget();
    int serveMatrixPropertyWidget();

    QStandardItemModel* tableModel;
    void setupProperties();
    void tableBtnPressed(int row);

    DataTypeWidget* basicPropertyWidget;
    MatrixPropertyWidget* matrixPropertyWidget;
    void matrixPropertyDetail(MaterialProperty* p);
};

#endif // MATERIALDIALOG_H
