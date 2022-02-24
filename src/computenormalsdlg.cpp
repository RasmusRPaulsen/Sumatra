#include "computenormalsdlg.h"
#include "ui_computenormalsdlg.h"

#include <sstream>

ComputeNormalsDlg::ComputeNormalsDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ComputeNormalsDlg)
{
    ui->setupUi(this);
}
//
//void ComputeNormalsDlg::Setup()
//{
//    ui->selectedSurface->addItem("item 1");
//    ui->selectedSurface->addItem("item 2");
//    ui->selectedSurface->addItem("item 3");
//    ui->selectedSurface->addItem("item 4");
//}

void ComputeNormalsDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
//    UpdateAllSceneData();
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

int ComputeNormalsDlg::getSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

bool ComputeNormalsDlg::replaceSource()
{
    return ui->ReplaceSourceChk->isChecked();
}

ComputeNormalsDlg::~ComputeNormalsDlg()
{
    delete ui;
}
