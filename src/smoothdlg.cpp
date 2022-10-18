#include "smoothdlg.h"
#include "ui_smoothdlg.h"

#include <sstream>

SmoothDlg::SmoothDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SmoothDlg)
{
    ui->setupUi(this);
}

SmoothDlg::~SmoothDlg()
{
    delete ui;
}


void SmoothDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void SmoothDlg::PopulateObjectCB()
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

int SmoothDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

bool SmoothDlg::GetReplaceSource()
{
    return ui->ReplaceSourceChk->isChecked();
}
