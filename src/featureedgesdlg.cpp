#include "featureedgesdlg.h"
#include "ui_featureedgesdlg.h"

#include <sstream>

FeatureEdgesDlg::FeatureEdgesDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FeatureEdgesDlg)
{
    ui->setupUi(this);
}

FeatureEdgesDlg::~FeatureEdgesDlg()
{
    delete ui;
}

void FeatureEdgesDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void FeatureEdgesDlg::PopulateObjectCB()
{
    if (!m3DScene)
        return;

    ui->selectedObjectCB->clear();
    for (int i = 0; i < m3DScene->GetNumberOfSurfaces(); i++)
    {
        std::ostringstream ost;
        ost << std::setfill('0') << std::setw(3) << i + 1 << " : " << m3DScene->GetSurfaceShortName(i);
        ui->selectedObjectCB->addItem(ost.str().c_str());
    }
    ui->selectedObjectCB->setCurrentIndex(0);
}

int FeatureEdgesDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

void FeatureEdgesDlg::getValues(bool& boundary, bool& nonManifold, bool& manifold, bool& sharp, double& sharpAngle)
{
    boundary = ui->boundaryEdgeChk->isChecked();
    nonManifold = ui->nonmanifoldChk->isChecked();
    manifold = ui->manifoldChk->isChecked();
    sharp = ui->sharpedgesChk->isChecked();
    sharpAngle = ui->sharpAngleSpin->value();
}
