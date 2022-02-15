#include "objectpropertiesdlg.h"
#include "ui_objectpropertiesdlg.h"

#include <sstream>

ObjectPropertiesDlg::ObjectPropertiesDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ObjectPropertiesDlg)
{
    ui->setupUi(this);

	connect(ui->selectedObjectCB, SIGNAL(currentIndexChanged(int)), SLOT(selectionChanged(int)));
}

ObjectPropertiesDlg::~ObjectPropertiesDlg()
{
    delete ui;
}

void ObjectPropertiesDlg::selectionChanged(int idx)
{
	UpdateAllSceneData();
}

void ObjectPropertiesDlg::Set3DScene(C3DScene* sc)
{
    m3DScene = sc;
	PopulateObjectCB();
	UpdateAllSceneData();
}

void ObjectPropertiesDlg::PopulateObjectCB()
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

void ObjectPropertiesDlg::UpdateAllSceneData()
{

	int idx = ui->selectedObjectCB->currentIndex();
	std::string fullName = m3DScene->GetSurfaceFullName(idx).c_str();
	ui->fullObjectName->setText(fullName.c_str());
	//ui->fullObjectName->
}
