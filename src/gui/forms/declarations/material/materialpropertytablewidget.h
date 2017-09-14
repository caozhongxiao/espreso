#ifndef MATERIALPROPERTYTABLEWIDGET_H
#define MATERIALPROPERTYTABLEWIDGET_H

#include <QWidget>
#include <QVector>
#include "../../../../config/configuration.h"
#include "../../elements/isavableobject.h"
#include "../datatypeeditwidget.h"

namespace espreso
{

    namespace Ui {
    class MaterialPropertyTableWidget;
    }

    class MaterialPropertyTableWidget : public QWidget, public ISavableObject
    {
        Q_OBJECT

    public:
        explicit MaterialPropertyTableWidget(QWidget *parent = 0, bool withHeader = true);
        ~MaterialPropertyTableWidget();

        void addProperty(ECFParameter* property);
        void addRow(const QString& name, ECFParameter* data, const QString& unit, const QString& symbol);

        void save() override;

    private:
        Ui::MaterialPropertyTableWidget *ui;

        QVector<DataTypeEditWidget*> m_rowWidgets;

        void createHeader();
    };

}
#endif // MATERIALPROPERTYTABLEWIDGET_H
