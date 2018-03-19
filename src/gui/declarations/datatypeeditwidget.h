#ifndef DATATYPEEDITWIDGET_H
#define DATATYPEEDITWIDGET_H

#include <QWidget>

#include <QLineEdit>
#include <QComboBox>
#include "tabletypewidget.h"
#include "piecewisetypewidget.h"
#include "../../config/configuration.h"
#include "../elements/expressionedit.h"
#include "../elements/ivalidatableobject.h"
#include "../elements/isavableobject.h"
#include "../elements/textitemwidget.h"

namespace espreso
{

    namespace Ui {
    class DataTypeEditWidget;
    }

    class DataTypeEditWidget : public TextItemWidget, public IValidatableObject,
            public ISavableObject
    {
        Q_OBJECT

    public:
        static QStringList typeNames();

        explicit DataTypeEditWidget(ECFParameter* data, QWidget* parent = 0);
        explicit DataTypeEditWidget(const std::vector<std::string>& variables, QWidget* parent = 0);
        ~DataTypeEditWidget();

        QComboBox* createComboBox(QWidget* parent = nullptr);
        void setComboBox(bool show);

        bool isValid() override;
        QString errorMessage() override;
        void save() override;

        virtual void setText(const QString& text) override;
        virtual QString text() override;

        QString value();
        void setValue(const QString& value);

        void setSharedDatatype(int *datatype);
        int datatype();

    private slots:
        void changeType(int index);

    private:
        DataTypeEditWidget(QWidget *parent = 0);

        Ui::DataTypeEditWidget *ui;

        ECFParameter* m_param;

        QComboBox* m_cmb;

        ExpressionEdit* uiExpression;
        TableTypeWidget* uiTable;
        PiecewiseTypeWidget* uiPiecewise;

        int activeType;
        int *m_shared = nullptr;
        std::vector<std::string> m_param_variables;
        QString param_getValue();
        std::vector<std::string> param_variables();

        void createUi();
        void initExpression();
        void initTable();
        void initPiecewise();
    };

}

#endif // DATATYPEEDITWIDGET_H
