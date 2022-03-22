#include "createprimitivedlg.h"
#include "ui_createprimitivedlg.h"

CreatePrimitiveDlg::CreatePrimitiveDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CreatePrimitiveDlg)
{
    ui->setupUi(this);
}

CreatePrimitiveDlg::~CreatePrimitiveDlg()
{
    delete ui;
}

void CreatePrimitiveDlg::setPickedPointPosition(double* p)
{
    ui->centerXSpn->setValue(p[0]);
    ui->centerYSpn->setValue(p[1]);
    ui->centerZSpn->setValue(p[2]);
}

int CreatePrimitiveDlg::getPrimitiveType()
{
    if (ui->sphereBtn->isChecked())
        return 0;
    else if (ui->cubeBtn->isChecked())
        return 1;
    return 0; // Sphere as default
}

void CreatePrimitiveDlg::getSphereParameters(double& r, double& phiStart, double& phiEnd, double& thetaStart, double& thetaEnd, 
    int& phiRes, int& thetaRes, double *center)
{
    r = ui->sphereRadiusSpn->value();
    phiStart = ui->spherePhiStartSpn->value();
    phiEnd = ui->spherePhiEndSpn->value();
    thetaStart = ui->sphereThetaStartSpn->value();
    thetaEnd = ui->sphereThetaEndSpn->value();
    phiRes = ui->spherePhiResSpn->value();
    thetaRes = ui->sphereThetaEndSpn->value();
    center[0] = ui->centerXSpn->value();
    center[1] = ui->centerYSpn->value();
    center[2] = ui->centerZSpn->value();
}
