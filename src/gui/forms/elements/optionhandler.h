#ifndef OPTIONHANDLER_H
#define OPTIONHANDLER_H

#include <QWidget>
#include "../../../config/configuration.h"

namespace espreso
{

    class OptionHandler : public QWidget
    {
        Q_OBJECT

    public:
        explicit OptionHandler(ECFParameter*, QWidget*);

    signals:
        void optionChanged();

    private slots:
        void onIndexChanged(int index);

    private:
        ECFParameter* m_option;
        bool optionsAdded = false;

    };

}

#endif // OPTIONHANDLER_H
