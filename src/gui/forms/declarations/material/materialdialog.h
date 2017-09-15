#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include <QDialog>
#include <QFrame>
#include <QVBoxLayout>
#include <QVector>
#include "../../../../config/configuration.h"
#include "../../elements/isavableobject.h"
#include "../../elements/ivalidatableobject.h"

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

        void accept() override;

    private slots:
        void redraw();

    private:
        Ui::MaterialDialog *ui;

        ECFObject* m_material;

        QFrame* m_frame;
        QVBoxLayout* m_frameLayout;

        QVector<ISavableObject*> m_save;
        QVector<IValidatableObject*> m_valid;

        void drawMe();
		void iterateObject(ECFObject*, QWidget* = nullptr);
        void drawOption(ECFParameter*, QWidget*);
    };

}

#endif // MATERIALDIALOG_H
