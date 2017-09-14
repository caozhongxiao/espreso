#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include <QDialog>
#include <QFrame>
#include <QVBoxLayout>
#include "../../../../config/configuration.h"

namespace espreso
{

    namespace Ui {
    class MaterialDialog;
    }

    class MaterialDialog : public QDialog
    {
        Q_OBJECT

    public:
        explicit MaterialDialog(ECFObject* material, QWidget *parent = 0);
        ~MaterialDialog();

    private slots:
        void redraw();

    private:
        Ui::MaterialDialog *ui;

        ECFObject* m_material;

        QFrame* m_frame;
        QVBoxLayout* m_frameLayout;

        void drawMe();
        void iterateObject(ECFObject*);
        void drawOption(ECFParameter*, QWidget*);
    };

}

#endif // MATERIALDIALOG_H
