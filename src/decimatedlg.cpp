#include "decimatedlg.h"
#include "ui_decimatedlg.h"

#include <sstream>


DecimateDlg::DecimateDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DecimateDlg)
{
    ui->setupUi(this);
}

DecimateDlg::~DecimateDlg()
{
    delete ui;
}


void DecimateDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void DecimateDlg::PopulateObjectCB()
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

int DecimateDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

bool DecimateDlg::GetReplaceSource()
{
    return ui->ReplaceSourceChk->isChecked();
}

