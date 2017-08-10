#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include <QDialog>
#include <QStandardItemModel>
#include "../data/material.h"
#include "../data/variable.h"
#include "datatypewidget.h"

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

    void visit(const BasicProperty& property) override;
    void visit(const IsotropicProperty& property) override;
    void visit(const DiagonalProperty& property) override;
    void visit(const SymmetricProperty& property) override;
    void visit(const AnisotropicProperty& property) override;

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
    void serveBasicPropertyWidget();
    void serveMatrixPropertyWidget();

    QStandardItemModel* tableModel;
    void setupProperties();
    void tableBtnPressed(int row);

    DataTypeWidget* basicPropertyWidget;
};

#endif // MATERIALDIALOG_H
