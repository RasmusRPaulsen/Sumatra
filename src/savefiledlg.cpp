#include "savefiledlg.h"
#include "ui_savefiledlg.h"

#include <sstream>

SaveFileDlg::SaveFileDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SaveFileDlg)
{
    ui->setupUi(this);
}

SaveFileDlg::~SaveFileDlg()
{
    delete ui;
}

void SaveFileDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void SaveFileDlg::PopulateObjectCB()
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

int SaveFileDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

bool SaveFileDlg::getApplyActorTransform()
{
    return ui->applyActorTransChk->isChecked();
}

bool SaveFileDlg::getWriteScalars()
{
    return ui->writeScalarsChk->isChecked();
}

bool SaveFileDlg::getWriteNormals()
{
    return ui->writeNormalsChk->isChecked();
}

bool SaveFileDlg::getWriteAsAscii()
{
    return ui->writeAsAsciiChk->isChecked();
}
