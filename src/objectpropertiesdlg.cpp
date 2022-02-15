#include "objectpropertiesdlg.h"
#include "ui_objectpropertiesdlg.h"

ObjectPropertiesDlg::ObjectPropertiesDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ObjectPropertiesDlg)
{
    ui->setupUi(this);
}

ObjectPropertiesDlg::~ObjectPropertiesDlg()
{
    delete ui;
}
