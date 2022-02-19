#include "objectpropertiesdlg.h"
#include "ui_objectpropertiesdlg.h"

#include <sstream>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <qcolordialog.h>



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
	connect(ui->RenderPointsAsSpheresChk, SIGNAL(toggled(bool)), SLOT(OnRenderingPointsAsSpheres(bool)));
	connect(ui->RenderLinesAsTubesChk, SIGNAL(toggled(bool)), SLOT(OnRenderingLinesAsTubes(bool)));
	connect(ui->BackFaceCullingChk, SIGNAL(toggled(bool)), SLOT(OnRenderingBackFaceCulling(bool)));
	connect(ui->FrontFaceCullingChk, SIGNAL(toggled(bool)), SLOT(OnRenderingFrontFaceCulling(bool)));
	connect(ui->EdgesVisibleChk, SIGNAL(toggled(bool)), SLOT(OnRenderingVisibleEdges(bool)));
	connect(ui->ShowScalarsChk, SIGNAL(toggled(bool)), SLOT(OnShowScalarsChk(bool)));
	connect(ui->PointSizeSpin, SIGNAL(valueChanged(int)), SLOT(OnPointSizeSpin(int)));
	connect(ui->LineWidthSpin, SIGNAL(valueChanged(int)), SLOT(OnLineWidthSpin(int)));
	connect(ui->OpacitySpin, SIGNAL(valueChanged(double)), SLOT(OnOpacitySpin(double)));
	connect(ui->AmbientSpin, SIGNAL(valueChanged(double)), SLOT(OnAmbientSpin(double)));
	connect(ui->DiffuseSpin, SIGNAL(valueChanged(double)), SLOT(OnDiffuseSpin(double)));
	connect(ui->SpecularSpin, SIGNAL(valueChanged(double)), SLOT(OnSpecularSpin(double)));
	connect(ui->SpecularPowerSpin, SIGNAL(valueChanged(double)), SLOT(OnSpecularPowerSpin(double)));
	connect(ui->ScalarRangeMin, SIGNAL(valueChanged(double)), SLOT(OnSetMinScalarRange(double)));
	connect(ui->ScalarRangeMax, SIGNAL(valueChanged(double)), SLOT(OnSetMaxScalarRange(double)));
	connect(ui->SpecularPowerSpin, SIGNAL(valueChanged(double)), SLOT(OnSpecularPowerSpin(double)));
	connect(ui->ChooseColorBtn, SIGNAL(clicked()), SLOT(OnChooseColorBtn()));
	connect(ui->DeleteObjectBtn, SIGNAL(clicked()), SLOT(OnDeleteObjectBtn()));
	connect(ui->RemoveScalarsBtn, SIGNAL(clicked()), SLOT(OnRemoveScalarsBtn()));
	connect(ui->RemoveNormalsBtn, SIGNAL(clicked()), SLOT(OnRemoveNormalsBtn()));
	connect(ui->ResetTransformBtn, SIGNAL(clicked()), SLOT(OnResetTransformBtn()));
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
	if (ui->selectedObjectCB->count() <= 0)
		return;

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
	ui->RenderPointsAsSpheresChk->setChecked(actor->GetProperty()->GetRenderPointsAsSpheres());
	ui->RenderLinesAsTubesChk->setChecked(actor->GetProperty()->GetRenderLinesAsTubes());
	ui->BackFaceCullingChk->setChecked(actor->GetProperty()->GetBackfaceCulling());
	ui->FrontFaceCullingChk->setChecked(actor->GetProperty()->GetFrontfaceCulling());
	ui->EdgesVisibleChk->setChecked(actor->GetProperty()->GetVertexVisibility());
	ui->PointSizeSpin->setValue(actor->GetProperty()->GetPointSize());
	ui->LineWidthSpin->setValue(actor->GetProperty()->GetLineWidth());
	ui->OpacitySpin->setValue(actor->GetProperty()->GetOpacity());
	ui->AmbientSpin->setValue(actor->GetProperty()->GetAmbient());
	ui->DiffuseSpin->setValue(actor->GetProperty()->GetDiffuse());
	ui->SpecularSpin->setValue(actor->GetProperty()->GetSpecular());
	ui->SpecularPowerSpin->setValue(actor->GetProperty()->GetSpecularPower());

	bool viewScalars = m3DScene->GetMapper(idx)->GetScalarVisibility();
	ui->ShowScalarsChk->setChecked(viewScalars);

	double minScals = 0;
	double maxScals = 0;
	if (viewScalars)
	{
		double range[2];
		m3DScene->GetScalarRange(idx, range);
		minScals = range[0];
		maxScals = range[1];
	}
	//ui->ScalarRangeMin->setText(QString::number(minScals));
	//ui->ScalarRangeMax->setText(QString::number(maxScals));
	ui->ScalarRangeMin->setValue(minScals);
	ui->ScalarRangeMax->setValue(maxScals);

	double color[3];
	actor->GetProperty()->GetColor(color);
	std::string info = m3DScene->GetSurfaceValues(idx);
	info = info + "Color: (" + std::to_string(color[0]) + ", " + std::to_string(color[1]) + ", " + std::to_string(color[2]) + ")\r\n";
	ui->ObjectInfo->setPlainText(info.c_str());

	// QPixmap pix = ui->ObjectColorLabel->pixmap();
	//int w = ui->ObjectColorLabel->rect().width();
	//int h = ui->ObjectColorLabel->rect().height();

	QPixmap pix(20, 20);
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

void ObjectPropertiesDlg::OnRenderingPointsAsSpheres(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetRenderPointsAsSpheres(on);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingLinesAsTubes(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetRenderLinesAsTubes(on);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingFrontFaceCulling(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetFrontfaceCulling(on);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingBackFaceCulling(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetBackfaceCulling(on);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRenderingVisibleEdges(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetEdgeVisibility(on);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnPointSizeSpin(int val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetPointSize(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnLineWidthSpin(int val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetLineWidth(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnOpacitySpin(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetOpacity(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnAmbientSpin(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetAmbient(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnDiffuseSpin(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetDiffuse(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnSpecularSpin(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetSpecular(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnSpecularPowerSpin(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->GetProperty()->SetSpecularPower(val);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnShowScalarsChk(bool on)
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetMapper(idx)->SetScalarVisibility(on);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnChooseColorBtn()
{
	int idx = ui->selectedObjectCB->currentIndex();
	double ocolor[3];
	m3DScene->GetActor(idx)->GetProperty()->GetColor(ocolor);

	QColor color = QColorDialog::getColor(QColor(ocolor[0] * 255, ocolor[1] * 255, ocolor[2] * 255), this);
	if (color.isValid())
	{
		m3DScene->GetActor(idx)->GetProperty()->SetColor(color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0);
		UpdateAllSceneData();
		emit valueChanged();
	}
}

// TODO: Do something when last object is deleted
void ObjectPropertiesDlg::OnDeleteObjectBtn()
{
	int idx = ui->selectedObjectCB->currentIndex();
	if (idx < 0)
		return;
	m3DScene->RemoveSurface(idx);
	PopulateObjectCB();
	UpdateAllSceneData();
	emit valueChanged();
}

void ObjectPropertiesDlg::OnResetTransformBtn()
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->GetActor(idx)->SetPosition(0, 0, 0);
	m3DScene->GetActor(idx)->SetScale(1);
	m3DScene->GetActor(idx)->SetOrientation(0, 0, 0);
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRemoveScalarsBtn()
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->RemoveScalars(idx);
	UpdateAllSceneData();
	emit valueChanged();
}

void ObjectPropertiesDlg::OnRemoveNormalsBtn()
{
	int idx = ui->selectedObjectCB->currentIndex();
	m3DScene->RemoveNormals(idx);
	UpdateAllSceneData();
	emit valueChanged();
}

void ObjectPropertiesDlg::OnSetMinScalarRange(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	double range[2];
	m3DScene->GetScalarRange(idx, range);
	if (val < range[1])
	{
		range[0] = val;
		m3DScene->SetScalarRange(idx, range);
		emit valueChanged();
	}
}

void ObjectPropertiesDlg::OnSetMaxScalarRange(double val)
{
	int idx = ui->selectedObjectCB->currentIndex();
	double range[2];
	m3DScene->GetScalarRange(idx, range);
	if (val > range[0])
	{
		range[1] = val;
		m3DScene->SetScalarRange(idx, range);
		emit valueChanged();
	}
}
