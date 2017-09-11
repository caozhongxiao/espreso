#include "formwidget.h"

#include <QLabel>
#include <QLineEdit>

using namespace espreso;

FormWidget::FormWidget(QWidget* parent) : QWidget(parent)
{
	this->m_layout = new QFormLayout;
	this->setLayout(m_layout);
}

void FormWidget::appendString(ECFParameter* p_string)
{
	QLabel* lbl = new QLabel(QString::fromStdString(p_string->metadata.description.at(0)),
							 this);

	QLineEdit* edit = new QLineEdit(QString::fromStdString(p_string->getValue()),
									this);

	this->m_layout->addRow(lbl, edit);
}
