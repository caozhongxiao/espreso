#ifndef FORMWIDGET_H
#define FORMWIDGET_H

#include <QWidget>
#include <QFormLayout>
#include "../../../config/configuration.h"

namespace espreso
{

	class FormWidget : public QWidget
	{
		Q_OBJECT

	public:
		explicit FormWidget(QWidget* = nullptr);

		void appendString(ECFParameter*);

	private:
		QFormLayout* m_layout;

	};

}

#endif // FORMWIDGET_H
