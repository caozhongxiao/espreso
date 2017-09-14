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

    class FormWidget : public QWidget, public ISavableObject
	{
		Q_OBJECT

	public:
		explicit FormWidget(QWidget* = nullptr);

		void appendString(ECFParameter*);

        void save() override;

	private:
		QFormLayout* m_layout;

        QVector<QPair<ECFParameter*, QLineEdit*> > m_strings;
	};

}

#endif // FORMWIDGET_H
