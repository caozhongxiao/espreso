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

        bool isValid() override;
        QString errorMessage() override;

	private:
		QFormLayout* m_layout;

        QVector<QPair<ECFParameter*, QLineEdit*> > m_strings;
        int m_invalidIndex = 0;
	};

}

#endif // FORMWIDGET_H
