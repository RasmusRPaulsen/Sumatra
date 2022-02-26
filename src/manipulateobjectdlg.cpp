#include "manipulateobjectdlg.h"
#include "ui_manipulateobjectdlg.h"

#include <sstream>

ManipulateObjectDlg::ManipulateObjectDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ManipulateObjectDlg)
{
    ui->setupUi(this);
}

ManipulateObjectDlg::~ManipulateObjectDlg()
{
    delete ui;
}

void ManipulateObjectDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
    PopulateObjectCB();
}

void ManipulateObjectDlg::PopulateObjectCB()
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

int ManipulateObjectDlg::GetSelectedSurface()
{
    return ui->selectedObjectCB->currentIndex();
}
