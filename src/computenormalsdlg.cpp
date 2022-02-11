#include "computenormalsdlg.h"
#include "ui_computenormalsdlg.h"

ComputeNormalsDlg::ComputeNormalsDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ComputeNormalsDlg)
{
    ui->setupUi(this);
}

void ComputeNormalsDlg::Setup()
{
    ui->selectedSurface->addItem("item 1");
    ui->selectedSurface->addItem("item 2");
    ui->selectedSurface->addItem("item 3");
    ui->selectedSurface->addItem("item 4");
}

int ComputeNormalsDlg::getSelectedSurface()
{
    return ui->selectedSurface->currentIndex();
}

bool ComputeNormalsDlg::replaceSource()
{
    return ui->replaceSource->isChecked();
}

ComputeNormalsDlg::~ComputeNormalsDlg()
{
    delete ui;
}
