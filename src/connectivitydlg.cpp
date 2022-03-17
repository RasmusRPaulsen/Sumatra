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
    setPickedPointPosition(sc->getLastPickedPoint());
}

void ConnectivityDlg::setPickedPointPosition(double* p)
{
    ui->pointXSpn->setValue(p[0]);
    ui->pointYSpn->setValue(p[1]);
    ui->pointZSpn->setValue(p[2]);
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

void ConnectivityDlg::getValues(bool& replaceSource, int& regionType, double* scalarRange, bool& fullScalarMode, double* point)
{
    replaceSource = ui->replaceSourceChk->isChecked();
    if (ui->allRegionsBtn->isChecked())
        regionType = 0;
    else if (ui->largestRegionBtn->isChecked())
        regionType = 1;
    else if (ui->outmostRegionBtn->isChecked())
        regionType = 2;
    else if (ui->closestPointRegionBtn->isChecked())
        regionType = 3;
    else if (ui->scalarConnectivityBtn->isChecked())
        regionType = 4;
    else
        regionType = 0; // TODO: Issue warning

    scalarRange[0] = ui->scalarMinSpn->value();
    scalarRange[1] = ui->scalarMaxSpn->value();
    fullScalarMode = ui->fullScalarConnectivityBtn->isChecked();
    point[0] = ui->pointXSpn->value();
    point[1] = ui->pointYSpn->value();
    point[2] = ui->pointZSpn->value();
}
