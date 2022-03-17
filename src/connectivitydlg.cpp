#include "connectivitydlg.h"
#include "ui_connectivitydlg.h"

#include <sstream>

ConnectivityDlg::ConnectivityDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ConnectivityDlg)
{
    ui->setupUi(this);
}

ConnectivityDlg::~ConnectivityDlg()
{
    delete ui;
}

void ConnectivityDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void ConnectivityDlg::PopulateObjectCB()
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

int ConnectivityDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}
