#include "formwidget.h"

#include <QLabel>

#include "../validators/validatorfactory.h"

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

    FieldHandler* edit = new FieldHandler(p_string);

    QPair<ECFParameter*, FieldHandler*> pair(p_string, edit);
    this->m_fields.append(pair);

    this->m_layout->addRow(lbl, edit);
}

void FormWidget::appendNonnegativeInteger(ECFParameter* p_nint)
{

    QLabel* lbl = new QLabel(QString::fromStdString(p_nint->metadata.description.at(0)),
                             this);

    ValidatorFactory* vf = new NonnegativeIntegerValidatorFactory;
    FieldHandler* edit = new FieldHandler(p_nint, vf);
    delete vf;

    QPair<ECFParameter*, FieldHandler*> pair(p_nint, edit);
    this->m_fields.append(pair);

    this->m_layout->addRow(lbl, edit);
}

void FormWidget::appendPositiveInteger(ECFParameter* p_nint)
{

    QLabel* lbl = new QLabel(QString::fromStdString(p_nint->metadata.description.at(0)),
                             this);

    ValidatorFactory* vf = new PositiveIntegerValidatorFactory;
    FieldHandler* edit = new FieldHandler(p_nint, vf);
    delete vf;

    QPair<ECFParameter*, FieldHandler*> pair(p_nint, edit);
    this->m_fields.append(pair);

    this->m_layout->addRow(lbl, edit);
}

void FormWidget::appendFloat(ECFParameter* p_nint)
{

    QLabel* lbl = new QLabel(QString::fromStdString(p_nint->metadata.description.at(0)),
                             this);

    ValidatorFactory* vf = new DoubleValidatorFactory;
    FieldHandler* edit = new FieldHandler(p_nint, vf);
    delete vf;

    QPair<ECFParameter*, FieldHandler*> pair(p_nint, edit);
    this->m_fields.append(pair);

    this->m_layout->addRow(lbl, edit);
}

void FormWidget::appendRow(const QString& label, QWidget* widget)
{
    this->m_layout->addRow(label, widget);
}

void FormWidget::save()
{
    for (auto pair = m_fields.cbegin();
         pair != m_fields.cend();
         ++pair)
    {
        pair->first->setValue(pair->second->value().toStdString());
    }
}

bool FormWidget::isValid()
{
    int index = 0;
    for (auto pair = m_fields.cbegin();
         pair != m_fields.cend();
         ++pair)
    {
        if (pair->second->value().isEmpty())
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
                     m_fields
                        .at(m_invalidIndex)
                        .first->metadata
                            .description
                            .at(0)
                     ));
}

void FormWidget::saveState()
{
    this->m_state_fields.clear();

    for (auto pair = m_fields.cbegin();
         pair != m_fields.cend();
         ++pair)
    {
        m_state_fields.append(
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
    for (auto pair = m_fields.cbegin();
         pair != m_fields.cend();
         ++pair)
    {
        pair->first->setValue(
                    m_state_fields.at(index)
                    );
        pair->second->setValue(
                    QString::fromStdString(m_state_fields.at(index))
                    );
        index++;
    }
}
