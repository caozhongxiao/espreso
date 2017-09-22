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

bool FormWidget::isValid()
{
    int index = 0;
    for (auto pair = m_strings.cbegin();
         pair != m_strings.cend();
         ++pair)
    {
        if (pair->second->text().isEmpty())
        {
            this->m_invalidIndex = index;
            return false;
        }
        index++;
    }
    return true;
}

QString FormWidget::errorMessage()
{
    return tr("Empty %1 field")
            .arg(QString::fromStdString(
                     m_strings
                        .at(m_invalidIndex)
                        .first->metadata
                            .description
                            .at(0)
                     ));
}

void FormWidget::saveState()
{
    this->m_state_strings.clear();

    for (auto pair = m_strings.cbegin();
         pair != m_strings.cend();
         ++pair)
    {
        m_state_strings.append(
                    pair->first->getValue()
                    );
    }

    this->m_stateStored = true;
}

void FormWidget::restoreState()
{
    if (!this->m_stateStored)
        return;

    int index = 0;
    for (auto pair = m_strings.cbegin();
         pair != m_strings.cend();
         ++pair)
    {
        pair->first->setValue(
                    m_state_strings.at(index)
                    );
        index++;
    }
}
