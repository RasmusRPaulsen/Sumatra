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
