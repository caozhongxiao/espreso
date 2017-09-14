#include "formwidget.h"

#include <QLabel>

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

    QPair<ECFParameter*, QLineEdit*> pair(p_string, edit);
    this->m_strings.append(pair);

	this->m_layout->addRow(lbl, edit);
}

void FormWidget::save()
{
    for (auto pair = m_strings.cbegin();
         pair != m_strings.cend();
         ++pair)
    {
        pair->first->setValue(pair->second->text().toStdString());
    }
}
