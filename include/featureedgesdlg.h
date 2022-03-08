#ifndef FEATUREEDGESDLG_H
#define FEATUREEDGESDLG_H

#include <QDialog>
#include "3DScene.h"


namespace Ui {
class FeatureEdgesDlg;
}

class FeatureEdgesDlg : public QDialog
{
    Q_OBJECT

public:
    explicit FeatureEdgesDlg(QWidget *parent = nullptr);
    ~FeatureEdgesDlg();


    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

    void getValues(bool& boundary, bool& nonManifold, bool& manifold, bool& sharp, double& sharpAngle);

private:
    C3DScene* m3DScene = NULL;

    Ui::FeatureEdgesDlg *ui;
};

#endif // FEATUREEDGESDLG_H
