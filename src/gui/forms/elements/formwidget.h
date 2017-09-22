#ifndef FORMWIDGET_H
#define FORMWIDGET_H

#include <QWidget>
#include <QFormLayout>
#include <QLineEdit>
#include "../../../config/configuration.h"
#include "isavableobject.h"
#include "ivalidatableobject.h"

namespace espreso
{

    class FormWidget : public QWidget, public ISavableObject,
            public IValidatableObject
	{
		Q_OBJECT

	public:
		explicit FormWidget(QWidget* = nullptr);

		void appendString(ECFParameter*);

        void save() override;
        void saveState() override;
        void restoreState() override;

        bool isValid() override;
        QString errorMessage() override;

	private:
		QFormLayout* m_layout;

        QVector<QPair<ECFParameter*, QLineEdit*> > m_strings;
        int m_invalidIndex = 0;

        QVector<std::string> m_state_strings;
        bool m_stateStored = false;
	};

}

#endif // FORMWIDGET_H
