#include "objectpropertiesdlg.h"
#include "ui_objectpropertiesdlg.h"

#include <sstream>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>


ObjectPropertiesDlg::ObjectPropertiesDlg(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ObjectPropertiesDlg)
{
    ui->setupUi(this);

	connect(ui->selectedObjectCB, SIGNAL(currentIndexChanged(int)), SLOT(selectionChanged(int)));
	connect(ui->RenderingPointsBtn, SIGNAL(toggled(bool)), SLOT(OnRenderingTypePoints(bool)));
	connect(ui->RenderingHiddenBtn, SIGNAL(toggled(bool)), SLOT(OnRenderingTypeHidden(bool)));
	connect(ui->RenderingWFBtn, SIGNAL(toggled(bool)), SLOT(OnRenderingTypeWireframe(bool)));
	connect(ui->RenderingSurfaceBtn, SIGNAL(toggled(bool)), SLOT(OnRenderingTypeSurface(bool)));
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

void ObjectPropertiesDlg::FullUpdate()
{
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

	vtkActor* actor = m3DScene->GetActor(idx);
	if (!actor)
		return;

	int rep = actor->GetProperty()->GetRepresentation();
	if (actor->GetVisibility() == 0)
	{
		ui->RenderingHiddenBtn->setChecked(true);
	}
	else if (rep == 0)
	{
		ui->RenderingPointsBtn->setChecked(true);
	}
	else if (rep == 1)
	{
		ui->RenderingWFBtn->setChecked(true);
	}
	else if (rep == 2)
	{
		ui->RenderingSurfaceBtn->setChecked(true);
	}
	double color[3];
	actor->GetProperty()->GetColor(color);
	std::string info = m3DScene->GetSurfaceValues(idx);
	info = info + "Color: (" + std::to_string(color[0]) + ", " + std::to_string(color[1]) + ", " + std::to_string(color[2]) + ")\r\n";
	ui->ObjectInfo->setPlainText(info.c_str());

	// QPixmap pix = ui->ObjectColorLabel->pixmap();
	QPixmap pix(16, 16);
	pix.fill(QColor((int)(color[0] * 255), (int)(color[1] * 255), (int)(color[2] * 255)));
	ui->ObjectColorLabel->setPixmap(pix);
}

void ObjectPropertiesDlg::OnRenderingTypePoints(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	if (on)
	{
		m3DScene->GetActor(idx)->GetProperty()->SetRepresentationToPoints();
		m3DScene->GetActor(idx)->VisibilityOn();
	}
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingTypeWireframe(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	if (on)
	{
		m3DScene->GetActor(idx)->GetProperty()->SetRepresentationToWireframe();
		m3DScene->GetActor(idx)->VisibilityOn();
	}
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingTypeSurface(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	if (on)
	{
		m3DScene->GetActor(idx)->GetProperty()->SetRepresentationToSurface();
		m3DScene->GetActor(idx)->VisibilityOn();
	}
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingTypeHidden(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	if (on)
	{
		m3DScene->GetActor(idx)->VisibilityOff();
	}
	emit valueChanged();
}

//
//void ObjectPropertiesDlg::OnRenderingTypeChange()
//{
//	//int idx = ui->selectedObjectCB->currentIndex();
//
//	//m3DScene->GetActor(idx)->SetVisibility(ui->RenderingHiddenBtn->isChecked());
//	//
//	//if (ui->RenderingPointsBtn->isChecked())
//	//	m3DScene->GetActor(idx)->GetProperty()->SetRepresentationToPoints();
//	//else if (ui->RenderingWFBtn->isChecked())
//	//	m3DScene->GetActor(idx)->GetProperty()->SetRepresentationToWireframe();
//	//else if (ui->RenderingSurfaceBtn->isChecked())
//	//	m3DScene->GetActor(idx)->GetProperty()->SetRepresentationToSurface();
//
//	//m3DScene->GetRenderer()->Render();
//}