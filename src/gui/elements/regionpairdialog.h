#ifndef REGIONPAIRDIALOG_H
#define REGIONPAIRDIALOG_H

#include <QDialog>

#include "../../config/configuration.h"
#include "../../config/ecf/ecf.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/plane/planeelement.h"

namespace espreso
{

namespace Ui {
class RegionPairDialog;
}

class RegionPairDialog : public QDialog
{
    Q_OBJECT

public:
    static RegionPairDialog* createRegionMaterial(ECFObject* map, Mesh* mesh, ECFObject* materials);
    static RegionPairDialog* createRegionExpression(ECFObject* map, Mesh* mesh);
    ~RegionPairDialog();

    void accept() override;

    QString region();

private:
    explicit RegionPairDialog(ECFDataType value,
                             ECFObject* map, Mesh* mesh,
                             ECFObject* scope, QWidget *parent = 0);

    Ui::RegionPairDialog *ui;

    Mesh* m_mesh;
    ECFObject* m_scope;
    ECFObject* m_map;
    ECFDataType m_first;
    ECFDataType m_second;
    QWidget* m_first_widget;
    QWidget* m_second_widget;

    QString m_region;

    QWidget* uiValue(ECFDataType type, QLayout* layout);

    std::string getKey();
};

}

#endif // REGIONPAIRDIALOG_H
