#include "computenormalsdlg.h"
#include "ui_computenormalsdlg.h"

#include <sstream>

ComputeNormalsDlg::ComputeNormalsDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ComputeNormalsDlg)
{
    ui->setupUi(this);
}

void ComputeNormalsDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void ComputeNormalsDlg::PopulateObjectCB()
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

int ComputeNormalsDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

bool ComputeNormalsDlg::GetReplaceSource()
{
    return ui->ReplaceSourceChk->isChecked();
}

bool ComputeNormalsDlg::GetFlipNormals()
{
    return ui->FlipNormalsChk->isChecked();
}

bool ComputeNormalsDlg::GetSplitNormals()
{
    return ui->SplitSharpEdgesChk->isChecked();
}

double ComputeNormalsDlg::GetSplitEdgeAngle()
{
    return ui->SharpEdgeAngleSpn->value();
}


ComputeNormalsDlg::~ComputeNormalsDlg()
{
    delete ui;
}
