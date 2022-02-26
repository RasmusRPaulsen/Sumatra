#include "manipulateobjectdlg.h"
#include "ui_manipulateobjectdlg.h"

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
