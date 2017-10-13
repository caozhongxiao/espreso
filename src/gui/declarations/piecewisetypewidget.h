#ifndef PIECEWISETYPEWIDGET_H
#define PIECEWISETYPEWIDGET_H

#include "tablewidget.h"
#include "../data/datatype.h"

namespace espreso
{

    class PiecewiseTypeWidget : public TableWidget
    {
        Q_OBJECT

    public:
        static QStringList headlines();

        PiecewiseTypeWidget(const std::vector<std::string>& variables,
                            QWidget* parent = 0);
        bool isValid() override;
        QString errorMessage() override;
        QString data() override;

        virtual void addData(const QString&) override;

    protected:
        virtual QString columnDefaultValue(int column) const override;

    private:
        std::vector<std::string> m_variables;
        QVector<QString> defaultValues;
        int m_invalidRow = 0;
        bool m_isValid = true;
    };

}

#endif // PIECEWISETYPEWIDGET_H
