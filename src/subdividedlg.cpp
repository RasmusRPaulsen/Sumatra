#include "subdividedlg.h"
#include "ui_subdividedlg.h"

#include <sstream>

SubdivideDlg::SubdivideDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SubdivideDlg)
{
    ui->setupUi(this);
}

SubdivideDlg::~SubdivideDlg()
{
    delete ui;
}


void SubdivideDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void SubdivideDlg::PopulateObjectCB()
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

int SubdivideDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}

bool SubdivideDlg::GetReplaceSource()
{
    return ui->ReplaceSourceChk->isChecked();
}

double SubdivideDlg::GetNumberOfSubdivisions()
{
    return ui->numSubDivSpn->value();
}

int SubdivideDlg::GetSubdivisionType()
{
    if (ui->linearRadio->isChecked())
        return 0;
    if (ui->butterflyRadio->isChecked())
        return 1;
    if (ui->loopRadio->isChecked())
        return 2;

    return 0;
}

