#include "3DScene.h"

//#include "CompatibilityFunctions.h"
#include "GeneralUtils.h"
#include "Process3DData.h"
// #include "RRPTestingColors.h"
#include "vtkPolyDataRandomResamplePoint.h"
//#include "vtkARANZReader.h"
//#include "vtkARANZWriter.h"
#include "vtkExtMisc.h"
#include "vtkPolyDataDifference.h"
#include "vtkRemoveUnusedPolyDataPoints.h"
#include "vtkTesselateBoundaryLoops.h"
#include "vtkTrapezoidSource.h"
#include "vtkPolyDataProjection.h"
//#include "MRFSurfaceReconstruction.h"
#include "GELRemeshing.h"
#include "vtkPolyDataPCANormals.h"
#include "vtkPolyDataOrientNormalsByVoting.h"
#include <vtkCamera.h>
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkAssembly.h>
#include <vtkAxesActor.h>
#include <vtkCaptionActor2D.h>
#include <vtkCellArray.h>
#include <vtkClipPolyData.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkDoubleArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLookUpTable.h>
#include <vtkMath.h>
#include <vtkPlane.h>
#include <vtkPlaneWidget.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPoints.h>
#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkProperty.h>
#include <vtkproperty2D.h>
#include <vtkTextProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkSphereSource.h>
#include <vtkSphereWidget.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkThresholdPoints.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkQuad.h>
#include <vtkPolyDataResamplePoints.h>
#include "CustomMultiTextureOBJImporter.h"

// #include <qmessagebox.h>

#include <sys/stat.h>
#include <algorithm>
#include <sstream>

C3DScene::C3DScene(CSumatraSettings* Settings)
{
	m_PointPicker = NULL;
	m_SphereMarkerLocator = NULL;
	m_Renderer = NULL;
	m_scalarBar = NULL;
	m_lookup = NULL;
	m_PlaneWidget = NULL;
	m_ScalarBarVisible = false;
	m_Axes = NULL;
	m_AxesVisible = false;
	m_MarkerValue = 1;
	m_StatusText = NULL;
	m_BackgroundActor = NULL;
	for (int i = 0; i < 6; i++)
		m_AllBounds[i] = 0;

	m_PickSphere = NULL;
	m_PickSphereMapper = NULL;
	m_PickSphereActor = NULL;
	m_MirrorID = -1;
	mSettings = Settings;
}

C3DScene::~C3DScene()
{
	if (m_Renderer)
		m_Renderer->Delete();

	if (m_scalarBar)
		m_scalarBar->Delete();

	if (m_lookup)
		m_lookup->Delete();

	if (m_PlaneWidget)
		m_PlaneWidget->Delete();

	if (m_Axes)
		m_Axes->Delete();

	if (m_SphereMarkerLocator)
		m_SphereMarkerLocator->Delete();

	if (m_StatusText)
		m_StatusText->Delete();

	if (m_BackgroundActor)
		m_BackgroundActor->Delete();

	for (unsigned int i = 0; i < m_Surfaces.size(); i++)
	{
		delete m_Surfaces[i];
	}

	if (m_PointPicker)
		m_PointPicker->Delete();

	if (m_PickSphere)
		m_PickSphere->Delete();
	if (m_PickSphereMapper)
		m_PickSphereMapper->Delete();
	if (m_PickSphereActor)
		m_PickSphereActor->Delete();

}

void C3DScene::Init()
{
	SetupScreenThings();
}

/*
# An example of how to set a gradient as the background of a render window.
# Goodwin Lawlor, 2006

set version [package require vtk]
wm withdraw .


# The background gradient Quad
vtkPoints quadPoints
quadPoints SetNumberOfPoints 4
quadPoints InsertPoint 0 0 0 0
quadPoints InsertPoint 1 1 0 0
quadPoints InsertPoint 2 1 1 0
quadPoints InsertPoint 3 0 1 0

vtkQuad aQuad
[aQuad GetPointIds] SetId 0 0
[aQuad GetPointIds] SetId 1 1
[aQuad GetPointIds] SetId 2 2
[aQuad GetPointIds] SetId 3 3

vtkUnsignedCharArray uchars
uchars SetNumberOfComponents 4
uchars SetNumberOfTuples 4
# grey to white
uchars SetTuple4 0 128 128 128 255;# bottom left RGBA colour
uchars SetTuple4 1 128 128 128 255;# bottom right RGBA colour
uchars SetTuple4 2 255 255 255 255;# top right RGBA colour
uchars SetTuple4 3 255 255 255 255;# top left RGBA colour
# blue to white
#   uchars SetTuple4 0 0 0 255 255;# bottom left RGBA colour
#   uchars SetTuple4 1 0 0 255 255;# bottom right RGBA colour
#   uchars SetTuple4 2 255 255 255 255;# top right RGBA colour
#   uchars SetTuple4 3 255 255 255 255;# top left RGBA colour
#try these by un-commenting them... 
#   uchars SetTuple4 0 255 0 0 255;# bottom left RGBA colour
#   uchars SetTuple4 1 0 0 255 255;# bottom right RGBA colour
#   uchars SetTuple4 2 255 0 0 255;# top right RGBA colour
#   uchars SetTuple4 3 0 255 0 255;# top left RGBA colour


vtkPolyData data
data Allocate 1 1
data InsertNextCell [aQuad GetCellType] [aQuad GetPointIds]
data SetPoints quadPoints
[data GetPointData] SetScalars uchars

vtkCoordinate coord
coord SetCoordinateSystemToNormalizedDisplay

vtkPolyDataMapper2D mapper2d
mapper2d SetInput data
mapper2d SetTransformCoordinate coord

vtkActor2D actor2d
actor2d SetMapper mapper2d
[actor2d GetProperty] SetDisplayLocationToBackground 


# a 3D object
vtkConeSource cone

vtkPolyDataMapper coneMapper
coneMapper SetInput [cone GetOutput]

vtkActor coneActor
coneActor SetMapper coneMapper

# Render the actors

vtkRenderer ren1
ren1 AddActor actor2d
ren1 AddActor coneActor

vtkRenderWindow renWin
renWin AddRenderer ren1


vtkInteractorStyleTrackballCamera style

vtkRenderWindowInteractor iren
iren SetRenderWindow renWin
iren SetInteractorStyle style
iren Initialize    

#if a (wish) console is available use it... otherwise use the vtk interactor
if {[string equal [info commands console] console]} {
iren AddObserver UserEvent {console show}
} else {
package require vtkinteraction
iren AddObserver UserEvent {wm deiconify .vtkInteract}
}


*/

void C3DScene::SetupScreenThings()
{
	m_Renderer = vtkRenderer::New();
	//m_Renderer->SetBackground(0.5, 0.5, 0.5);
	m_Renderer->SetBackground(mSettings->mBackgroundColor);

	m_lookup = vtkLookupTable::New();
	SetScalarLookupTableNum(0);

	m_Axes = vtkAxesActor::New();
	m_Axes->SetVisibility(m_AxesVisible);
	m_Axes->SetConeResolution(100);
//	m_Axes->SetShaftTypeToCylinder();
//	m_Renderer->AddActor(m_Axes);

	m_StatusText = vtkTextActor::New();
	m_StatusText->SetInput("Hello world!");
	m_StatusText->GetTextProperty()->SetFontSize(18);
	m_StatusText->GetTextProperty()->SetFontFamilyToArial();
	m_StatusText->VisibilityOff();

		//SetFontSize 18
		//$tprop SetFontFamilyToArial
		//$tprop SetJustificationToCentered


	m_Renderer->AddActor2D(m_StatusText);

	// Gradient background
	// From http://www.bioengineering-research.com/vtk/BackgroundGradient.tcl
	vtkPoints *quadPoints = vtkPoints::New();
	quadPoints->SetNumberOfPoints(4);
	quadPoints->InsertPoint(0,0,0,0);
	quadPoints->InsertPoint(1,1,0,0);
	quadPoints->InsertPoint(2,1,1,0);
	quadPoints->InsertPoint(3,0,1,0);

	vtkQuad *aQuad = vtkQuad::New();
	aQuad->GetPointIds()->SetId(0, 0);
	aQuad->GetPointIds()->SetId(1, 1);
	aQuad->GetPointIds()->SetId(2, 2);
	aQuad->GetPointIds()->SetId(3, 3);

	vtkUnsignedCharArray *uchars = vtkUnsignedCharArray::New();
	uchars->SetNumberOfComponents(4);
	uchars->SetNumberOfTuples(4);
	uchars->SetTuple4(0, 128, 128, 128, 255); // bottom left RGBA colour
	uchars->SetTuple4(1, 128, 128, 128, 255); // bottom right RGBA colour
	uchars->SetTuple4(2, 255, 255, 255, 255); // top right RGBA colour
	uchars->SetTuple4(3, 255, 255, 255, 255); // top left RGBA colour

//# blue to white
//#   uchars SetTuple4 0 0 0 255 255;# bottom left RGBA colour
//#   uchars SetTuple4 1 0 0 255 255;# bottom right RGBA colour
//#   uchars SetTuple4 2 255 255 255 255;# top right RGBA colour
//#   uchars SetTuple4 3 255 255 255 255;# top left RGBA colour
//#try these by un-commenting them... 
//#   uchars SetTuple4 0 255 0 0 255;# bottom left RGBA colour
//#   uchars SetTuple4 1 0 0 255 255;# bottom right RGBA colour
//#   uchars SetTuple4 2 255 0 0 255;# top right RGBA colour
//#   uchars SetTuple4 3 0 255 0 255;# top left RGBA colour

	vtkPolyData *bgPD = vtkPolyData::New();
	bgPD->Allocate(1, 1);
	bgPD->InsertNextCell(aQuad->GetCellType(), aQuad->GetPointIds());
	bgPD->SetPoints(quadPoints);
	bgPD->GetPointData()->SetScalars(uchars);

	vtkCoordinate *bgCoord = vtkCoordinate::New();
	bgCoord->SetCoordinateSystemToNormalizedDisplay();

	vtkPolyDataMapper2D *bgMapper = vtkPolyDataMapper2D::New();
	bgMapper->SetInputData(bgPD);
	bgMapper->SetTransformCoordinate(bgCoord);

	m_BackgroundActor = vtkActor2D::New();
	m_BackgroundActor->SetMapper(bgMapper);
	m_BackgroundActor->GetProperty()->SetDisplayLocationToBackground();

	m_Renderer->AddActor(m_BackgroundActor);

	quadPoints->Delete();
	aQuad->Delete();
	uchars->Delete();
	bgPD->Delete();
	bgCoord->Delete();
	bgMapper->Delete();

	m_BackgroundActor->SetVisibility(false);

	SetupScalarBar();

	m_PointPicker = vtkPointPicker::New();
	m_PointPicker->SetTolerance(0.01);

	m_PickSphere = vtkSphereSource::New();
	m_PickSphere->SetThetaResolution(50);
	m_PickSphere->SetPhiResolution(50);

	m_PickSphereMapper = vtkPolyDataMapper::New();
	m_PickSphereMapper->SetInputConnection(m_PickSphere->GetOutputPort());

	m_PickSphereActor = vtkActor::New();
	m_PickSphereActor->SetMapper(m_PickSphereMapper);
	m_PickSphereActor->SetPickable(0);
	m_PickSphereActor->SetVisibility(0);

	m_Renderer->AddActor(m_PickSphereActor);
}

void C3DScene::RemoveAllSurfaces()
{
	for (unsigned int i = 0; i < m_Surfaces.size(); i++)
	{
		m_Renderer->RemoveActor(m_Surfaces[i]->m_actor);
		delete m_Surfaces[i];
	}
	m_Surfaces.clear();
	UpdateAllBounds();
}

void C3DScene::RemoveSurface(unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, QString("Warning"), QString("Selected object not valid!"));
		return;
	}

	m_Renderer->RemoveActor(m_Surfaces[objID]->m_actor);
	
	CSurfaceProperties *tp =  m_Surfaces[objID];

	std::vector<CSurfaceProperties*>::iterator ti = m_Surfaces.begin() + objID;
//		&m_Surfaces[objID];
	m_Surfaces.erase(ti);

	delete tp;

	UpdateAllBounds();
}


void C3DScene::InitialisePlaneWidget(vtkRenderWindowInteractor *interactor, unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	m_PlaneWidget = vtkPlaneWidget::New();
	m_PlaneWidget->SetResolution(50);
	m_PlaneWidget->SetRepresentationToSurface();
	m_PlaneWidget->SetInteractor(interactor);
	m_PlaneWidget->GetPlaneProperty()->SetOpacity(0.8);
	m_PlaneWidget->GetPlaneProperty()->SetColor(0.0, 0.8, 0.8);
	m_PlaneWidget->SetInputData(m_Surfaces[objID]->m_polyData);
    m_PlaneWidget->SetPlaceFactor(1.00);
    m_PlaneWidget->PlaceWidget();
	m_PlaneWidget->On();
}

void C3DScene::InitialiseSphereWidget(vtkRenderWindowInteractor *interactor, unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	m_SphereWidget = vtkSphereWidget::New();
	m_SphereWidget->SetPhiResolution(50);
	m_SphereWidget->SetThetaResolution(50);

	m_SphereWidget->SetRepresentationToSurface();
	m_SphereWidget->SetInteractor(interactor);
	m_SphereWidget->GetSphereProperty()->SetOpacity(0.8);
	m_SphereWidget->GetSphereProperty()->SetColor(0.0, 0.8, 0.8);
	m_SphereWidget->SetInputData(m_Surfaces[objID]->m_polyData);
	m_SphereWidget->SetPlaceFactor(0.15);
	m_SphereWidget->ScaleOff();
	m_SphereWidget->PlaceWidget();
	m_SphereWidget->On();
}

void C3DScene::DeleteSphereWidget()
{
	if (m_SphereWidget)
	{
		m_SphereWidget->Off();
		m_SphereWidget->SetInteractor(NULL);
		m_SphereWidget->Delete();

		m_SphereWidget = NULL;
	}
}

void C3DScene::DeletePlaneWidget()
{
	if (m_PlaneWidget)
	{
		m_PlaneWidget->Off();
		m_PlaneWidget->SetInteractor(NULL);
		m_PlaneWidget->Delete();

		m_PlaneWidget = NULL;
	}
}

void C3DScene::ResetCamera()
{
	m_Renderer->ResetCamera();
}

vtkRenderer *C3DScene::GetRenderer() const
{
	return m_Renderer;
}


void C3DScene::AddSurfaceToRenderer(CSurfaceProperties *sp)
{
	m_Surfaces.push_back(sp);
	m_Renderer->AddActor(sp->m_actor);
	if (sp->HaveScalars())
	{
		if (sp->m_polyData->GetPointData()->GetScalars()->GetNumberOfComponents() == 1)
			SetScalarBarVisible(true);
	}
	UpdateAllBounds();
}

bool C3DScene::ReadFile(const std::string &fname)
{
	//! Very special hack for multi-surface obj files
	if (CGeneralUtils::GetExtensionFromFilename(fname) == "obj")
	{
		CCustomMultiTextureOBJImporter OBJImporter;
		if (OBJImporter.ReadFile(fname))
		{
			for (int i = 0; i < OBJImporter.NumberOfMeshes(); i++)
			{
				std::ostringstream ost;
				ost <<  CGeneralUtils::StripPathAndExtensionFromFilename(fname) << "_"  << i;
				CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);

				surfProbs->InitialiseSurface(OBJImporter.GetMesh(i));
				if (OBJImporter.GetTexture(i) != NULL)
					surfProbs->SetTexture(OBJImporter.GetTexture(i));
				surfProbs->m_shortname = ost.str();

				AddSurfaceToRenderer(surfProbs);
			}

			// Only reset camera first time
			if (m_Surfaces.size() == OBJImporter.NumberOfMeshes())
			{
				ResetCamera();
			}

			return true;
		}
	}
	
	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);

	if (!surfProbs->ReadFromFile(fname))
	{
		delete surfProbs;
		return false;
	}

	AddSurfaceToRenderer(surfProbs);

	// Only reset camera first time
	if (m_Surfaces.size() == 1)
	{
		ResetCamera();
	}

	return true;
}

bool C3DScene::SaveFile(unsigned int id, const std::string &fname, bool ApplyUserTransform, bool Binary, bool WriteNormals, bool WriteScalars)
{
	if (id >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return false;
	}

	return m_Surfaces[id]->SaveToFile(fname, ApplyUserTransform, Binary, WriteNormals, WriteScalars);
}


int C3DScene::GetNumberOfSurfaces() const
{
	return m_Surfaces.size();
}

void C3DScene::FlipPlane()
{
	if (m_PlaneWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No plane!");
		return;
	}

	double n[3];
	m_PlaneWidget->GetNormal(n);
	m_PlaneWidget->SetNormal(-n[0], -n[1], -n[2]);
}

void C3DScene::MirrorWithPlane( unsigned int objID )
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	if (m_PlaneWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No plane!");
		return;
	}

	// Use actor transform before cutting
	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(m_Surfaces[objID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(m_Surfaces[objID]->m_polyData);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	std::string mirrorName = m_Surfaces[objID]->m_shortname + "_mirror";

	bool MirrorExist = false;
	// Check if mirrored version already exists (by checking the name)
	if (m_MirrorID >= 0 && m_MirrorID < (int)m_Surfaces.size())
	{
		if (m_Surfaces[m_MirrorID]->m_shortname == mirrorName)
		{
			MirrorExist = true;
		}
	}

	vtkPolyData *mirror = vtkPolyData::New();
	mirror->DeepCopy(ptransSource->GetOutput());

	double O[3];
	double N[3];
	m_PlaneWidget->GetOrigin(O);
	m_PlaneWidget->GetNormal(N);

	vtkDataArray *normals = mirror->GetPointData()->GetNormals();


	for (int i = 0; i < mirror->GetNumberOfPoints(); i++)
	{
		double p[3];
		mirror->GetPoint(i, p);
	
		// Vector from origon to point
		double v[3];
		v[0] = p[0] - O[0];
		v[1] = p[1] - O[1];
		v[2] = p[2] - O[2];
		
		// Distance of point from plane
		double dist = vtkMath::Dot(v, N);

		// New mirrored point
		double pp[3];
		pp[0] = p[0] - 2 * dist * N[0];
		pp[1] = p[1] - 2 * dist * N[1];
		pp[2] = p[2] - 2 * dist * N[2];
		
		mirror->GetPoints()->SetPoint(i, pp);

		if (normals)
		{
			double nn[3];
			normals->GetTuple(i, nn);

			// Transform tip of normal
			double ppp[3];
			ppp[0] = p[0] - nn[0];
			ppp[1] = p[1] - nn[1];
			ppp[2] = p[2] - nn[2];

			// Vector from origon to tip of normal
			v[0] = ppp[0] - O[0];
			v[1] = ppp[1] - O[1];
			v[2] = ppp[2] - O[2];

			// Distance of point from plane
			dist = vtkMath::Dot(v, N);

			// New mirrored tip of normal
			double pppp[3];
			pppp[0] = ppp[0] - 2 * dist * N[0];
			pppp[1] = ppp[1] - 2 * dist * N[1];
			pppp[2] = ppp[2] - 2 * dist * N[2];

			nn[0] = pppp[0]-pp[0];
			nn[1] = pppp[1]-pp[1];
			nn[2] = pppp[2]-pp[2];

			vtkMath::Normalize(nn);

			normals->SetTuple(i, nn);
		}
	}

	if (MirrorExist)
	{
		m_Surfaces[m_MirrorID]->m_polyData->DeepCopy(mirror);
	}
	else
	{
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs2->InitialiseSurface(mirror);
		surfProbs2->m_shortname = mirrorName;

		AddSurfaceToRenderer(surfProbs2);
		m_MirrorID = m_Surfaces.size() - 1;
	}

	ptransSource->Delete();
	sourceTrans->Delete();
}


void C3DScene::CutWithPlane(unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	if (m_PlaneWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No plane!");
		return;
	}

	// Use actor transform before cutting
	vtkTransform *sourceTrans = vtkTransform::New();
	 sourceTrans->SetMatrix(m_Surfaces[objID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	 ptransSource->SetInputData(m_Surfaces[objID]->m_polyData);
	 ptransSource->SetTransform(sourceTrans);
	 ptransSource->Update();
	
	vtkPlane *cplane = vtkPlane::New();
	cplane->SetOrigin(m_PlaneWidget->GetOrigin());
	cplane->SetNormal(m_PlaneWidget->GetNormal());

	vtkClipPolyData *clipper = vtkClipPolyData::New();
	clipper->SetInputConnection(ptransSource->GetOutputPort());
	clipper->SetClipFunction(cplane);
	clipper->GenerateClipScalarsOff();
	clipper->GenerateClippedOutputOff();
	clipper->SetValue(0);
	clipper->InsideOutOff();
	clipper->Update();

//	vtkPolyDataConnectivityFilter *connect = vtkPolyDataConnectivityFilter::New();
//	 connect->SetInput(clipper->GetOutput());
//	 connect->SetExtractionModeToLargestRegion();
//	 connect->Update();


	m_Surfaces[objID]->StoreState();

//	m_Surfaces[objID]->m_polyData->DeepCopy(connect->GetOutput());
	m_Surfaces[objID]->m_polyData->DeepCopy(clipper->GetOutput());
	
	m_Surfaces[objID]->m_actor->SetPosition(0,0,0);
	m_Surfaces[objID]->m_actor->SetScale(1);
	m_Surfaces[objID]->m_actor->SetOrientation(0,0,0);

	std::string cname = m_Surfaces[objID]->m_shortname;
	// This is so clumsy
	if (cname.size() < 3 || cname[cname.size()-3] != 'c' || cname[cname.size()-2] != 'u' || cname[cname.size()-1] != 't')
	{
		m_Surfaces[objID]->m_shortname += "_cut";
	}

//	connect->Delete();
	clipper->Delete();
	cplane->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
}

void C3DScene::UndoLastOperation(unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	m_Surfaces[objID]->RestoreState();
}


void C3DScene::StartMarkWithSphere(unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	if (m_SphereWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No Sphere!");
		return;
	}

	// Check if object has scalars. If not generate some
	if (m_Surfaces[objID]->m_polyData->GetPointData()->GetScalars() == NULL)
	{
		vtkDoubleArray *scalars = vtkDoubleArray::New();
		scalars->SetNumberOfComponents(1);
		scalars->SetNumberOfValues(m_Surfaces[objID]->m_polyData->GetNumberOfPoints());

		double val = 0;
		for (int i = 0; i < scalars->GetNumberOfTuples(); i++)
		{
			scalars->SetTuple(i, &val);
		}

		m_Surfaces[objID]->m_polyData->GetPointData()->SetScalars(scalars);
		scalars->Delete();

		m_Surfaces[objID]->UpdateScalarProperties();
	}

	// Use actor transform before marking
	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(m_Surfaces[objID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(m_Surfaces[objID]->m_polyData);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	m_SphereMarkerLocator = vtkPointLocator::New();
	m_SphereMarkerLocator->SetDataSet(ptransSource->GetOutput());
	m_SphereMarkerLocator->SetNumberOfPointsPerBucket(1);
	m_SphereMarkerLocator->BuildLocator();

	sourceTrans->Delete();
	ptransSource->Delete();

	m_StatusText->VisibilityOn();
}

void C3DScene::EndMarkWithSphere(unsigned int objID)
{
	if (m_SphereMarkerLocator)
	{
		m_SphereMarkerLocator->Delete();
		m_SphereMarkerLocator = NULL;
	}
	m_StatusText->VisibilityOff();
}

void C3DScene::MarkWithSphere(unsigned int objID)
{
	if (objID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	if (m_SphereWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No Sphere!");
		return;
	}

	double cp[3];
	m_SphereWidget->GetCenter(cp);
	double rad = m_SphereWidget->GetRadius();

	vtkIdList *pids = vtkIdList::New();

	m_SphereMarkerLocator->FindPointsWithinRadius (rad, cp, pids);

	vtkPolyData *pd = m_Surfaces[objID]->m_polyData;

	for (int i = 0; i < pids->GetNumberOfIds(); i++)
	{
		vtkIdType id = pids->GetId(i);

		pd->GetPointData()->GetScalars()->SetTuple(id, &m_MarkerValue);
	}
	pd->GetPointData()->GetScalars()->Modified();
	pd->Modified();
	
	if (pids->GetNumberOfIds() > 0)
	{
		double scale[2];
		GetScalarRange(objID, scale);
		if (m_MarkerValue < scale[0] || m_MarkerValue > scale[1])
		{
			if (m_MarkerValue < scale[0])
				scale[0] = m_MarkerValue;
			if (m_MarkerValue > scale[1])
				scale[1] = m_MarkerValue;
			SetScalarRange(objID, scale);
		}
		m_Surfaces[objID]->UpdateScalarProperties();
	}
	pids->Delete();
}


void C3DScene::FlipObject(unsigned int SourceID, bool ReplaceSurface)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	// Use actor transform before cutting
	vtkTransform *sourceTrans = vtkTransform::New();
	 sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	 ptransSource->SetInputData(m_Surfaces[SourceID]->m_polyData);
	 ptransSource->SetTransform(sourceTrans);
	 ptransSource->Update();

	vtkTransform *trans = vtkTransform::New();
	 trans->Scale(-1, 1, 1);

	vtkTransformPolyDataFilter *filt = vtkTransformPolyDataFilter::New();
	 filt->SetInputConnection(ptransSource->GetOutputPort());
	 filt->SetTransform(trans);
	 filt->Update();

 	vtkDataArray *normals = filt->GetOutput()->GetPointData()->GetNormals();

    bool hasNormals = (normals != NULL);

	if (ReplaceSurface)
	{
		if (hasNormals)
		{
			vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
			 norms->SetInputConnection(filt->GetOutputPort());
			 norms->SetSplitting(0);
			 norms->ConsistencyOn();
			 norms->SetFeatureAngle(90);
			 norms->Update();

			m_Surfaces[SourceID]->m_polyData->DeepCopy(norms->GetOutput());
			norms->Delete();
		}
		else
		{
			m_Surfaces[SourceID]->m_polyData->DeepCopy(filt->GetOutput());
		}
	}
	else
	{
		std::string name = m_Surfaces[SourceID]->m_shortname + "_flip";
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);

		if (hasNormals)
		{
			vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
			 norms->SetInputConnection(filt->GetOutputPort());
			 norms->SetSplitting(0);
			 norms->ConsistencyOn();
			 norms->SetFeatureAngle(90);
			 norms->Update();

			surfProbs2->InitialiseSurface(norms->GetOutput());
			norms->Delete();
		}
		else
		{
			surfProbs2->InitialiseSurface(filt->GetOutput());
		}
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}
	filt->Delete();
	trans->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
}


void C3DScene::VisualiseNormals(unsigned int SourceID, bool VisualiseAdjacentFlipped)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkDataArray *normals = m_Surfaces[SourceID]->m_polyData->GetPointData()->GetNormals();
	if (!normals)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No normals computed for object!");
		return;
	}

	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(m_Surfaces[SourceID]->m_polyData);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	normals = ptransSource->GetOutput()->GetPointData()->GetNormals();

	// Get an estimate of the size of the pointcloud
	double bounds[6];
	ptransSource->GetOutput()->GetBounds(bounds);

	double l = sqrt((bounds[1]-bounds[0]) * (bounds[1]-bounds[0]) + 
		(bounds[3]-bounds[2]) * (bounds[3]-bounds[2]) +
		(bounds[5]-bounds[4]) * (bounds[5]-bounds[4]));

	// Normal length. 1 % of bounding box diagonal
	double NL = l * 0.01;

	vtkPolyData *pd = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *lines = vtkCellArray::New();

	for (int i = 0; i < ptransSource->GetOutput()->GetNumberOfPoints(); i++)
	{
		double p[3];
		double n[3];
		double pt[3];

		ptransSource->GetOutput()->GetPoint(i, p);
		normals->GetTuple(i, n);

		vtkMath::Normalize(n);
		n[0] *= NL;
		n[1] *= NL;
		n[2] *= NL;

		pt[0] = n[0] + p[0];
		pt[1] = n[1] + p[1];
		pt[2] = n[2] + p[2];

		lines->InsertNextCell(2);
		vtkIdType id = points->InsertNextPoint(p);
		lines->InsertCellPoint(id);
		id = points->InsertNextPoint(pt);
		lines->InsertCellPoint(id);
	}

	pd->SetPoints(points);
	points->Delete();
	pd->SetLines(lines);
	lines->Delete();

	std::string name = m_Surfaces[SourceID]->m_shortname + "_LineNormals";
	CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);

	surfProbs2->InitialiseSurface(pd);
	surfProbs2->m_shortname = name;

	AddSurfaceToRenderer(surfProbs2);

	if (VisualiseAdjacentFlipped)
	{
		vtkPolyData *pd2 = vtkPolyData::New();
		int nFlipped = CGELRemeshing::VisualiseAdjacentFlippedNormals(ptransSource->GetOutput(), pd2, NL);
		if (nFlipped == 0)
		{
			// TODO: Update QMessageBox::information(this, "Warning","No adjacent faces with flipped normals");
		}
		else
		{ 
			name = m_Surfaces[SourceID]->m_shortname + "_AdjacentFlippedNormals";
			CSurfaceProperties *surfProbs3 = new CSurfaceProperties(m_lookup, mSettings);

			surfProbs3->InitialiseSurface(pd2);
			surfProbs3->m_shortname = name;

			AddSurfaceToRenderer(surfProbs3);
		}

		pd2->Delete();
	}


	pd->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
}

void C3DScene::CloseHoles(unsigned int SourceID, bool ReplaceSurface)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	// Use actor transform before cutting
	vtkTransform *sourceTrans = vtkTransform::New();
	 sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	 ptransSource->SetInputData(m_Surfaces[SourceID]->m_polyData);
	 ptransSource->SetTransform(sourceTrans);
	 ptransSource->Update();

	vtkTesselateBoundaryLoops *tess = vtkTesselateBoundaryLoops::New();
	 tess->SetInputConnection(ptransSource->GetOutputPort());
	 tess->AppendTesselationToInputOn();
	 tess->Update();

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(tess->GetOutput());
	}
	else
	{
		std::string name = m_Surfaces[SourceID]->m_shortname + "_CloseHoles";
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);

		surfProbs2->InitialiseSurface(tess->GetOutput());
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}
	tess->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
}

//
//void C3DScene::CreateFastRBFSurface(unsigned int SourceID, bool OutputSurface, bool OutputEstNormals, bool OutputDensity,
//							double MinNormalLength, double MaxNormalLength, double Accuracy, double Resolution,
//							double EstNormalsRadius, double EstNormalsPlaneFactor, bool ErrorBarFitting,
//							double ISOSmooth, bool DefineDensity, int InsideID, int OutsideID)
//{
//	if (SourceID >= m_Surfaces.size())
//	{
//		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
//		return;
//	}
//	
//	const std::string FastRBF = "\"C:\\Program Files\\FarField Technology\\FastRBF v1.4\\Command Line\\FastRBF.exe\"";
//	struct stat fs;
//	if (stat(FastRBF.c_str(), &fs) == 0)
//	{
//		// TODO: Update QMessageBox::information(this, "Warning","FastRBF not installed on this computer!");
//		return;
//	}
//
////	double accuracy   = 0.2;
////	double resolution = 1.0;
//
//	char *tempvar = getenv("TEMP");
//	
//	if( tempvar == NULL )
//	{
//		// TODO: Update QMessageBox::information(this, "Warning","TEMP variable not set");
//		return;
//	}
//	
//	std::string sname = m_Surfaces[SourceID]->m_shortname;
//
//	std::string outname = std::string(tempvar) + "\\" + sname;
//	std::string onameRaw = outname + ".aranz";
//	std::string onameNormals = outname + ".RBFEstimatedNormals.aranz";
//	std::string onameDensity = outname + ".Density.aranz";
//	std::string onameRBF = outname + ".RBF.aranz";
//	std::string onameISO = outname + ".ISOSurf.aranz";
//
//	std::ostringstream cmdline;
//
//	if (DefineDensity)
//	{
//		if (SourceID == InsideID || SourceID == OutsideID || InsideID == OutsideID)
//		{
//			// TODO: Update QMessageBox::information(this, "Warning","Three different surfaces should be selected");
//			return;
//		}
//		vtkPoints *pts = vtkPoints::New();
//		vtkCellArray *verts = vtkCellArray::New();
//		vtkPolyData *pd = vtkPolyData::New();
//		vtkDoubleArray *scalars = vtkDoubleArray::New();
//		scalars->SetNumberOfComponents(1);
//
//		int j;
//		for (j = 0; j < m_Surfaces[InsideID]->m_polyData->GetNumberOfPoints(); j++)
//		{
//			double val = -1;
//			double p[3];
//			m_Surfaces[InsideID]->m_polyData->GetPoint(j, p);
//
//			vtkIdType id = pts->InsertNextPoint(p);
//			verts->InsertNextCell(1);
//			verts->InsertCellPoint(id);
//			scalars->InsertTuple(id, &val);
//		}
//		for (j = 0; j < m_Surfaces[OutsideID]->m_polyData->GetNumberOfPoints(); j++)
//		{
//			double val = 1;
//			double p[3];
//			m_Surfaces[OutsideID]->m_polyData->GetPoint(j, p);
//
//			vtkIdType id = pts->InsertNextPoint(p);
//			verts->InsertNextCell(1);
//			verts->InsertCellPoint(id);
//			scalars->InsertTuple(id, &val);
//		}
//		for (j = 0; j < m_Surfaces[SourceID]->m_polyData->GetNumberOfPoints(); j++)
//		{
//			double val = 0;
//			double p[3];
//			m_Surfaces[SourceID]->m_polyData->GetPoint(j, p);
//
//			vtkIdType id = pts->InsertNextPoint(p);
//			verts->InsertNextCell(1);
//			verts->InsertCellPoint(id);
//			scalars->InsertTuple(id, &val);
//		}
//
//		pd->SetPoints(pts);
//		pd->SetVerts(verts);
//		pd->GetPointData()->SetScalars(scalars);
//
//		vtkARANZWriter *writer = vtkARANZWriter::New();
//		writer->SetInputData(pd);
//		writer->SetFileName(onameDensity.c_str());
//		writer->Write();
//		writer->Delete();
//
//		scalars->Delete();
//		verts->Delete();
//		pts->Delete();
//		pd->Delete();
//	}
//
//	if (!DefineDensity)
//	{
//		// Convert file to something that aranz software can read
//		vtkARANZWriter *writer = vtkARANZWriter::New();
//		writer->SetInputData(m_Surfaces[SourceID]->m_polyData);
//		writer->SetFileName(onameRaw.c_str());
//		writer->Write();
//		writer->Delete();
//	
//		// Estimate normals
//		cmdline << FastRBF << " estimatenormals -ascii -radius=" << EstNormalsRadius 
//			<< " -factor=" << EstNormalsPlaneFactor << " -points " << onameRaw << " " << onameNormals;
//		system(cmdline.str().c_str());
//		cmdline.str("");
//		
//		// Density from normals
//		cmdline << FastRBF << " densityfromnormals -ascii -length=" << MinNormalLength << "," << MaxNormalLength << " " << onameNormals << " " << onameDensity;
//		system(cmdline.str().c_str());
//		cmdline.str("");
//	}
//
//	// Unique
//	cmdline << FastRBF << " unique -ascii " << onameDensity << " " << onameDensity;
//	system(cmdline.str().c_str());
//	cmdline.str("");
//
//	// Fitting RBF
//	if (!ErrorBarFitting)
//		cmdline << FastRBF << " fit -reduce -accuracy=" << Accuracy << " " << onameDensity << " " << onameRBF;
//	else
//		cmdline << FastRBF << " fit -errorbar -accuracy=" << Accuracy << " " << onameDensity << " " << onameRBF;
//
//	std::cout << "Trying:\n" << cmdline.str() << std::endl;
//	system(cmdline.str().c_str());
//	cmdline.str("");
//
//	// ISO surfacing
//	if (ISOSmooth > 0)
//	{
//		cmdline << FastRBF << " isosurf -ascii -resolution=" << Resolution << " -smooth=" << ISOSmooth << " " << onameRBF << " " << onameISO;
//	}
//	else
//	{
//		cmdline << FastRBF << " isosurf -ascii -resolution=" << Resolution << " " << onameRBF << " " << onameISO;
//	}
//	std::cout << "Trying:\n" << cmdline.str() << std::endl;
//	system(cmdline.str().c_str());
//	cmdline.str("");
//
//	if (OutputSurface)
//	{
//		vtkARANZReader *ISOSurf = vtkARANZReader::New();
//		 ISOSurf->SetFileName(onameISO.c_str());
//		 ISOSurf->Update();
//		if (ISOSurf->GetOutput()->GetNumberOfPoints() > 0)
//		{
//			CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
//			surfProbs->InitialiseSurface(ISOSurf->GetOutput());
//			surfProbs->m_shortname = m_Surfaces[SourceID]->m_shortname + "_FastRBF";
//
//			AddSurfaceToRenderer(surfProbs);
//		}
//		ISOSurf->Delete();
//	}
//	if (OutputEstNormals && !DefineDensity)
//	{
//		vtkARANZReader *reader = vtkARANZReader::New();
//		reader->SetFileName(onameNormals.c_str());
//		reader->Update();
//
//		if (reader->GetOutput()->GetNumberOfPoints() > 0)
//		{
//			CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
//			surfProbs->InitialiseSurface(reader->GetOutput());
//			surfProbs->m_shortname = m_Surfaces[SourceID]->m_shortname + "_FastRBFEstimatedNormals";
//
//			AddSurfaceToRenderer(surfProbs);
//		}
//
//		reader->Delete();
//	}
//	if (OutputDensity)
//	{
//		vtkARANZReader *reader = vtkARANZReader::New();
//		reader->SetFileName(onameDensity.c_str());
//		reader->Update();
//
//		if (reader->GetOutput()->GetNumberOfPoints() > 0)
//		{
//			CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
//			surfProbs->InitialiseSurface(reader->GetOutput());
//			surfProbs->m_shortname = m_Surfaces[SourceID]->m_shortname + "_FastRBFDensity";
//
//			AddSurfaceToRenderer(surfProbs);
//		}
//
//		reader->Delete();
//	}
//
//	// Clean up
//	remove(onameRaw.c_str());
//	remove(onameNormals.c_str());
//	remove(onameDensity.c_str());
//	remove(onameRBF.c_str());
//	remove(onameISO.c_str());
//}

void C3DScene::SmoothSurface(unsigned int SourceID, bool ReplaceSurface, int smoothtype, int NumIt, double RelaxFactor,
			bool BoundarySmooth, bool FeatureEdgeSmooth, double FeatureAngle, bool GenerateErrScal)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkPolyData *smooth = vtkPolyData::New();
	CProcess3DData::DoSmoothSurface(m_Surfaces[SourceID]->m_polyData, smooth, smoothtype, NumIt, RelaxFactor, BoundarySmooth, FeatureEdgeSmooth, FeatureAngle, GenerateErrScal);

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(smooth);
	}
	else
	{
		std::ostringstream name;
		name << m_Surfaces[SourceID]->m_shortname;
		if (smoothtype == 0)
			name << "_LaplacianSmooth";
		else
			name << "_SincSmooth";

		name << "_it_" << NumIt << "_relax_" << RelaxFactor << "_BoundarySmooth_" << BoundarySmooth
			<< "_FeatureEdgeSmooth_" << FeatureEdgeSmooth << "_FeatureAngle_" << FeatureAngle;

		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs->InitialiseSurface(smooth);
		surfProbs->m_shortname = name.str();

		AddSurfaceToRenderer(surfProbs);
	}

	smooth->Delete();
}


void C3DScene::DecimateSurface(unsigned int SourceID, bool ReplaceSurface, int decimtype, float decimfactor, bool preservetopology)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkPolyData *decim = vtkPolyData::New();
	CProcess3DData::DoDecimateSurface(m_Surfaces[SourceID]->m_polyData, decim, decimtype, decimfactor, preservetopology);

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(decim);
	}
	else
	{
		std::ostringstream name;
		name << m_Surfaces[SourceID]->m_shortname;
		if (decimtype == 0)
			name << "_DecimatePro_" << decimfactor;
		if (decimtype == 1)
			name << "_QuadricDecimation_" << decimfactor;

		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs->InitialiseSurface(decim);
		surfProbs->m_shortname = name.str();

		AddSurfaceToRenderer(surfProbs);
	}

	decim->Delete();
}

void C3DScene::ComputeFeatureEdges(unsigned int SourceID, int EdgeType, double FeatureAngle)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	vtkPolyData *FeatureEdges = vtkPolyData::New();
	CProcess3DData::DoComputeFeatureEdges(m_Surfaces[SourceID]->m_polyData, FeatureEdges, EdgeType, FeatureAngle);

	std::ostringstream name;
	name << m_Surfaces[SourceID]->m_shortname;
	if (EdgeType == 0)
		name << "_BoundaryEdges";
	if (EdgeType == 1)
		name << "_NonManifoldEdges";
	if (EdgeType == 2)
		name << "_ManifoldEdges";
	if (EdgeType == 3)
		name << "_FeatureEdges_" << FeatureAngle;

	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs->InitialiseSurface(FeatureEdges);
	surfProbs->m_shortname = name.str();

	AddSurfaceToRenderer(surfProbs);

	FeatureEdges->Delete();
}

void C3DScene::SubdivideSurface(unsigned int SourceID, bool ReplaceSurface, int subdivtype, int subdivisions)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	vtkPolyData *subdiv = vtkPolyData::New();
	CProcess3DData::DoSubdivideSurface(m_Surfaces[SourceID]->m_polyData, subdivisions, subdiv, subdivtype);

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(subdiv);
	}
	else
	{
		std::ostringstream name;
		name << m_Surfaces[SourceID]->m_shortname;
		if (subdivtype == 0)
			name << "_LinearSubdivided_" << subdivisions;
		if (subdivtype == 1)
			name << "_ButterflySubdivided_" << subdivisions;
		if (subdivtype == 2)
			name << "_LoopSubdivided_" << subdivisions;

		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs->InitialiseSurface(subdiv);
		surfProbs->m_shortname = name.str();

		AddSurfaceToRenderer(surfProbs);
	}

	subdiv->Delete();
}

void C3DScene::Connectivity(unsigned int SourceID, bool ReplaceSurface, int RegionType)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	if (RegionType == 0)
	{
		vtkPolyDataConnectivityFilter *connect = vtkPolyDataConnectivityFilter::New();
		connect->SetInputData(m_Surfaces[SourceID]->m_polyData);
		connect->ColorRegionsOff();
		connect->SetExtractionModeToSpecifiedRegions();
		connect->Update();

		int N = connect->GetNumberOfExtractedRegions();

		std::ostringstream nn;

		for (int i = 0; i < N; i++)
		{
			connect->InitializeSpecifiedRegionList();
			connect->AddSpecifiedRegion(i);
			connect->Update();
			vtkRemoveUnusedPolyDataPoints *rem = vtkRemoveUnusedPolyDataPoints::New();
			rem->SetInputConnection(connect->GetOutputPort());
			rem->Update();

			nn.str("");
			nn << m_Surfaces[SourceID]->m_shortname << "_connect" << i;
			CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
			surfProbs->InitialiseSurface(rem->GetOutput());
			surfProbs->m_shortname = nn.str();

			AddSurfaceToRenderer(surfProbs);
			rem->Delete();
		}
		connect->Delete();
	}
	else if (RegionType == 1)
	{
		vtkPolyDataConnectivityFilter *connect = vtkPolyDataConnectivityFilter::New();
		connect->SetInputData(m_Surfaces[SourceID]->m_polyData);
		connect->ColorRegionsOff();
		connect->SetExtractionModeToLargestRegion();
		connect->Update();

		vtkRemoveUnusedPolyDataPoints *rem = vtkRemoveUnusedPolyDataPoints::New();
		rem->SetInputConnection(connect->GetOutputPort());
		rem->Update();

		if (ReplaceSurface)
		{
			m_Surfaces[SourceID]->m_polyData->DeepCopy(rem->GetOutput());
		}
		else
		{
			std::string name = m_Surfaces[SourceID]->m_shortname + "_largestregion";
			CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
			surfProbs->InitialiseSurface(rem->GetOutput());
			surfProbs->m_shortname = name;
			AddSurfaceToRenderer(surfProbs);
		}
		connect->Delete();
		rem->Delete();
	}	
	else if (RegionType == 2)
	{
		vtkPolyData *con = vtkPolyData::New();
		CProcess3DData::ExtractOuterSurface(m_Surfaces[SourceID]->m_polyData, con);

		if (ReplaceSurface)
		{
			m_Surfaces[SourceID]->m_polyData->DeepCopy(con);
		}
		else
		{
			std::string name = m_Surfaces[SourceID]->m_shortname + "_outersurface";
			CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
			surfProbs->InitialiseSurface(con);
			surfProbs->m_shortname = name;
			AddSurfaceToRenderer(surfProbs);
		}
		con->Delete();
	}
}


void C3DScene::CalculateMRFSurface( unsigned int SourceID, bool ReplaceSurface, bool RecomputeNorms, int priortype, 
	bool LargestCliqueOnly, bool MarchingCubes, bool ConComp, int InputType, int NumberOfPoints, int PPerNormal, int PPerDistance, double TriangleSizeFactor, bool UseTargetEdgeLengths )
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
//
//
//	CMRFSurfaceReconstruction::CParms parms;
//	//parms.InputName = inputname;
//	//parms.InputDir = inputdir;
//	//parms.OutPutDir = outputdir;
//	//	parms.CalibrationDir = calibrationDir;
//	parms.SetInputType(InputType);
//	parms.EstimateNormals = RecomputeNorms;
//	parms.ConnectedComponentAnalysis = ConComp;
//	parms.LargestCliqueOnly = LargestCliqueOnly;
//	parms.WriteLevel = 0;
////	parms.SampleFactor = 0.10;
//	parms.NormalRadius = 2.5;
//	if (PPerNormal > 0)
//	{
//		parms.NormalSearchMode = 1;
//		parms.NormalPointsPerNormal	 = PPerNormal;
//	}
//	parms.NumberOfDistances = PPerDistance;
//	parms.DistanceMode = 2;
//	parms.Polygoniser = 0;
//	if (MarchingCubes)
//		parms.Polygoniser = 1;
//	parms.GlobalBeta = 0.9;
//	parms.WeightMode = 1;
//	parms.LocalWeightMaxDist = 3.0;
//	parms.PriorType = priortype;
//	parms.Optimisation = 0;
//	parms.MaxVolumeSize = 10000000;
////	parms.EnlargePercent = 200;
//	parms.PadVoxels = 5;
//
//	parms.MaxIterations = 1000;
//	parms.CellSizeFactor = TriangleSizeFactor;
//	parms.UseLocalWeights = true;
//	parms.BandedICM = true;
//	parms.AddNoise = false;
//	parms.NoiseNature = 1;
//	parms.UseLeavePatchOut = false;
//	parms.MultiLevel = true;
//	parms.Remesh = true;
//	parms.RemeshToTargetEdgeLengths = UseTargetEdgeLengths;
//	parms.AdaptiveParameters = true;
//	parms.AggresiveCrop = true;
//
//	std::string msg;
//	CMRFSurfaceReconstruction MRFSurface;
//	bool result = true;
//	
//	// check for normals
//	if (parms.InputType == 1 ||  parms.InputType ==  2 || parms.InputType == 4 || parms.InputType ==  5 )
//	{
//		vtkDataArray *normals = m_Surfaces[SourceID]->m_polyData->GetPointData()->GetNormals(); 
//		if (!normals)
//		{
//			// TODO: Update QMessageBox::information(this, "Warning","Selected input does not have normals");
//			return;
//		}
//	}
//
//	if (NumberOfPoints < 0)
//	{
//		// use all points
//		result = MRFSurface.OneShotCompute(parms, m_Surfaces[SourceID]->m_polyData, msg);
//	}
//	else
//	{
//		// reduce input point set size
//		vtkPolyDataResamplePoints *resampler = vtkPolyDataResamplePoints::New();
//		resampler->SetInputData(m_Surfaces[SourceID]->m_polyData);
//		resampler->SetNumberOfOutputPoints(NumberOfPoints);
//		resampler->Update();
//
//		result = MRFSurface.OneShotCompute(parms, resampler->GetOutput(), msg);
//		resampler->Delete();
//	}
//	if (!result)
//	{
//		// TODO: Update QMessageBox::information(this, "Warning",msg.c_str());
//		return;
//	}
//
//	if (ReplaceSurface)
//	{
//		m_Surfaces[SourceID]->m_polyData->DeepCopy(MRFSurface.GetMRFSurfaceCropped());
//		// Remove scalar information
//		m_Surfaces[SourceID]->m_polyData->GetPointData()->SetScalars(NULL); 
//	}
//	else
//	{
//		std::string name = m_Surfaces[SourceID]->m_shortname + "_MRFSurface";
//		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
//		surfProbs->InitialiseSurface(MRFSurface.GetMRFSurface());
//		surfProbs->m_shortname = name;
//		AddSurfaceToRenderer(surfProbs);
//		surfProbs->m_actor->SetVisibility(0);
//
//		name = m_Surfaces[SourceID]->m_shortname + "_MRFSurface_Cropped";
//		CSurfaceProperties *surfProbs1 = new CSurfaceProperties(m_lookup, mSettings);
//		surfProbs1->InitialiseSurface(MRFSurface.GetMRFSurfaceCropped());
//		surfProbs1->m_polyData->GetPointData()->SetScalars(NULL);
//		surfProbs1->m_shortname = name;
//
//		AddSurfaceToRenderer(surfProbs1);
//
//		name = m_Surfaces[SourceID]->m_shortname + "_MRFSurface_AggresiveCrop";
//		CSurfaceProperties *surfProbs3 = new CSurfaceProperties(m_lookup, mSettings);
//		surfProbs3->InitialiseSurface(MRFSurface.GetMRFSurfaceAgressivelyCropped());
//		surfProbs3->m_shortname = name;
//
//		AddSurfaceToRenderer(surfProbs3);
//		surfProbs3->m_actor->SetVisibility(0);
//
//		name = m_Surfaces[SourceID]->m_shortname + "_NormalData";
//		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
//		surfProbs2->InitialiseSurface(MRFSurface.GetNormalData());
//		surfProbs2->m_shortname = name;
//
//		AddSurfaceToRenderer(surfProbs2);
//		surfProbs2->m_actor->SetVisibility(0);
//	}
}

void C3DScene::CalculateNormals(unsigned int SourceID, bool ReplaceSurface, bool flip, bool split, double featureAngle)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
	 norms->SetInputData(m_Surfaces[SourceID]->m_polyData);
	 norms->SetSplitting(split);
	 norms->SetFeatureAngle(featureAngle);
	 norms->SetFlipNormals(flip);
	 norms->ConsistencyOn();
	 norms->NonManifoldTraversalOff();
	 norms->Update();

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(norms->GetOutput());
	}
	else
	{
		std::string name = m_Surfaces[SourceID]->m_shortname + "_normals";
		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs->InitialiseSurface(norms->GetOutput());
		surfProbs->m_shortname = name;

		AddSurfaceToRenderer(surfProbs);
	}
	norms->Delete();
}

void C3DScene::ScalarFiltering(unsigned int SourceID, bool ReplaceSurface, double SphereRadius)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	if (m_Surfaces[SourceID]->m_polyData->GetPointData()->GetScalars() == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object does not contain scalar values");
		return;
	}

	// Oldscalars
	vtkDoubleArray *oscals = vtkDoubleArray::SafeDownCast(m_Surfaces[SourceID]->m_polyData->GetPointData()->GetScalars());

	vtkPolyData *pd = vtkPolyData::New();
	pd->DeepCopy(m_Surfaces[SourceID]->m_polyData);

	// Trying to smooth  scalars
	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(pd);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	vtkDoubleArray *nscals = vtkDoubleArray::New();
	nscals->SetNumberOfTuples(oscals->GetNumberOfTuples());

	for (int i = 0; i < oscals->GetNumberOfTuples(); i++)
	{
		double p[3];
		pd->GetPoint(i, p);

		vtkIdList *neighPts = vtkIdList::New();

		locator->FindPointsWithinRadius(SphereRadius, p, neighPts);

		int nIds = neighPts->GetNumberOfIds();
		if (nIds != 0)
		{
			// Should always be true
			double scalSum = 0;
			for (int n = 0; n < nIds; n++)
			{
				vtkIdType cid = neighPts->GetId(n);
				scalSum += oscals->GetValue(cid);
			}
			scalSum /= (double)nIds;
			nscals->SetValue(i, scalSum);
		}
		else
		{
			// Should not happen
			nscals->SetValue(i, oscals->GetValue(i));
		}
		neighPts->Delete();
	}
	pd->GetPointData()->SetScalars(nscals);
	nscals->Delete();

	std::ostringstream ost;
	ost.str("");
	ost << m_Surfaces[SourceID]->m_shortname << "_ScalarSmooth_" << SphereRadius;

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(pd);
		m_Surfaces[SourceID]->m_shortname = ost.str();
		m_Surfaces[SourceID]->UpdateScalarProperties();
	}
	else
	{
		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);

		surfProbs->m_shortname = ost.str();

		vtkTransform *sourceTrans = vtkTransform::New();
		sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

		vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
		ptransSource->SetInputData(pd);
		ptransSource->SetTransform(sourceTrans);
		ptransSource->Update();

		surfProbs->InitialiseSurface(ptransSource->GetOutput());

		AddSurfaceToRenderer(surfProbs);

		ptransSource->Delete();
		sourceTrans->Delete();
	}

	pd->Delete();
	locator->Delete();
}

void C3DScene::ScalarThreshold(unsigned int SourceID, bool ReplaceSurface, double lowThreshold, double highThreshold)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	if (m_Surfaces[SourceID]->m_polyData->GetPointData()->GetScalars() == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object does not contain scalar values");
		return;
	}

 	vtkClipPolyData *clipper = vtkClipPolyData::New();
 	clipper->SetInputData(m_Surfaces[SourceID]->m_polyData);
 	clipper->GenerateClipScalarsOff();
 	clipper->GenerateClippedOutputOff();
 	clipper->SetValue(lowThreshold);
 	clipper->InsideOutOff();
 	clipper->Update();

	vtkClipPolyData *clipper2 = vtkClipPolyData::New();
 	clipper2->SetInputConnection(clipper->GetOutputPort());
 	clipper2->GenerateClipScalarsOff();
 	clipper2->GenerateClippedOutputOff();
 	clipper2->SetValue(highThreshold);
 	clipper2->InsideOutOn();
 	clipper2->Update();


	if(clipper2->GetOutput()->GetNumberOfPoints() == 0)
	{
	  // TODO: Update QMessageBox::information(this, "Warning","Thresholding resulted in zero points, new point cloud was not created!");
	  clipper->Delete();
	  clipper2->Delete();
	}
	else
	{
	  std::ostringstream ost;
	  ost.str("");
	  ost << m_Surfaces[SourceID]->m_shortname << "_threshold_" << lowThreshold << "_to_" << highThreshold;

      if (ReplaceSurface)
	  {
		m_Surfaces[SourceID]->m_polyData->DeepCopy(clipper2->GetOutput());
		m_Surfaces[SourceID]->m_shortname = ost.str();
		m_Surfaces[SourceID]->UpdateScalarProperties();
	  }
	  else
	  {
		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
		
		surfProbs->m_shortname = ost.str();

		vtkTransform *sourceTrans = vtkTransform::New();
		sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

		vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
		ptransSource->SetInputConnection(clipper2->GetOutputPort());
	    ptransSource->SetTransform(sourceTrans);
	    ptransSource->Update();
	
		surfProbs->InitialiseSurface(ptransSource->GetOutput());

		AddSurfaceToRenderer(surfProbs);

		ptransSource->Delete();
		sourceTrans->Delete();
	  }
	}
	clipper->Delete();
	clipper2->Delete();
}


void C3DScene::DoMerge(unsigned int SourceID, unsigned int TargetID)
{
	if (SourceID >= m_Surfaces.size() || 
		TargetID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected objects not valid!");
		return;
	}

	vtkPolyData *source = m_Surfaces[SourceID]->m_polyData;
	vtkPolyData *target = m_Surfaces[TargetID]->m_polyData;
	vtkMatrix4x4 *sourcePreTransform = m_Surfaces[SourceID]->m_actor->GetMatrix();
	vtkMatrix4x4 *targetPreTransform = m_Surfaces[TargetID]->m_actor->GetMatrix();

	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(sourcePreTransform);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(source);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	vtkTransform *targetTrans = vtkTransform::New();
	targetTrans->SetMatrix(targetPreTransform);

	vtkTransformPolyDataFilter *ptransTarget = vtkTransformPolyDataFilter::New();
	ptransTarget->SetInputData(target);
	ptransTarget->SetTransform(targetTrans);
	ptransTarget->Update();


	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkPolyData *pd = vtkPolyData::New();

	vtkDataArray *normals1 = ptransSource->GetOutput()->GetPointData()->GetNormals(); 
	vtkDataArray *normals2 = ptransTarget->GetOutput()->GetPointData()->GetNormals(); 

	vtkDoubleArray *Normals = NULL;
	if (normals1 && normals2)
	{
		Normals = vtkDoubleArray::New();
		Normals->SetNumberOfComponents(3);
		Normals->SetName("Normals");

	}

	for (int i = 0; i < ptransSource->GetOutput()->GetNumberOfPoints(); i++)
	{
		double p[3];
		ptransSource->GetOutput()->GetPoint(i, p);

		vtkIdType id = pts->InsertNextPoint(p);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);

		if (Normals)
		{
			double n[3];
			normals1->GetTuple(i, n);
			Normals->InsertNextTuple(n);
		}
	}

	for (int i = 0; i < ptransTarget->GetOutput()->GetNumberOfPoints(); i++)
	{
		double p[3];
		ptransTarget->GetOutput()->GetPoint(i, p);

		vtkIdType id = pts->InsertNextPoint(p);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);

		if (Normals)
		{
			double n[3];
			normals2->GetTuple(i, n);
			Normals->InsertNextTuple(n);
		}
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->SetNormals(Normals);
	pts->Delete();
	verts->Delete();

	if (Normals)
	{
		Normals->Delete();
	}

	std::string name = m_Surfaces[SourceID]->m_shortname + "_MergedWith_" + m_Surfaces[TargetID]->m_shortname;

	CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs2->InitialiseSurface(pd);
	surfProbs2->m_shortname = name;

	AddSurfaceToRenderer(surfProbs2);

	pd->Delete();
	ptransTarget->Delete();
	targetTrans->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
}



bool C3DScene::DoCompare( unsigned int SourceID, unsigned int TargetID, bool SignedDistance, bool ReplaceSource, bool ExcludeEdgePoints )
{
	if (SourceID >= m_Surfaces.size() || 
		TargetID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected objects not valid!");
		return false;
	}

	vtkPolyData *source = m_Surfaces[SourceID]->m_polyData;
	vtkPolyData *target = m_Surfaces[TargetID]->m_polyData;
	vtkMatrix4x4 *sourcePreTransform = m_Surfaces[SourceID]->m_actor->GetMatrix();
	vtkMatrix4x4 *targetPreTransform = m_Surfaces[TargetID]->m_actor->GetMatrix();

	if (target->GetNumberOfPolys() < 1)
	{
		// TODO: Update QMessageBox::information(this, "Warning","Target must be a triangulated surface");
		return false;
	}

	// Check for normals. If signed compare and no normals for the target, then compute them
	vtkDataArray *normals = target->GetPointData()->GetNormals(); 

	vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
	if (!normals)
	{
		norms->SetInputData(target);
		norms->SplittingOff();
		norms->ConsistencyOn();
		norms->Update();

		target = norms->GetOutput();
	}

	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(sourcePreTransform);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(source);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	vtkTransform *targetTrans = vtkTransform::New();
	targetTrans->SetMatrix(targetPreTransform);

	vtkTransformPolyDataFilter *ptransTarget = vtkTransformPolyDataFilter::New();
	ptransTarget->SetInputData(target);
	ptransTarget->SetTransform(targetTrans);
	ptransTarget->Update();

	vtkPolyDataDifference *PDDiff = vtkPolyDataDifference::New();
	PDDiff->SetSignedDistance(SignedDistance);
	PDDiff->SetExcludeEdgePoints(ExcludeEdgePoints);

	PDDiff->SetInputConnection(ptransSource->GetOutputPort());
	PDDiff->SetTargetData(ptransTarget->GetOutput());
	PDDiff->Update();

	std::string name = m_Surfaces[SourceID]->m_shortname + "_ComparedTo_"	+ m_Surfaces[TargetID]->m_shortname;

	if (ReplaceSource)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(PDDiff->GetOutput());
		m_Surfaces[SourceID]->m_shortname = name;
		m_Surfaces[SourceID]->UpdateScalarProperties();
		SetScalarBarVisible(true);
	}
	else
	{
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs2->InitialiseSurface(PDDiff->GetOutput());
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}
	PDDiff->Delete();
	ptransTarget->Delete();
	targetTrans->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
	norms->Delete();

	return true;
}

void C3DScene::DoRemesh(unsigned int SourceID, bool ReplaceSurface, bool ScalarTargetLengths, double UniformTargetLength)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	if (ScalarTargetLengths && m_Surfaces[SourceID]->m_polyData->GetPointData()->GetScalars() == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object does not contain scalar values");
		return;
	}

	if (ScalarTargetLengths)
	{
		vtkPolyData *pdLengths = vtkPolyData::New();
		pdLengths->DeepCopy(m_Surfaces[SourceID]->m_polyData);

		vtkPolyData *pd = vtkPolyData::New();

		CGELRemeshing remesher;

		bool res = remesher.RemeshWithTargetLengthsDirectlyOnMesh(m_Surfaces[SourceID]->m_polyData, pdLengths, pd);

		if (!res)
		{
			// TODO: Update QMessageBox::information(this, "Warning","Something went wrong in remeshing");
		}
		else
		{
			std::ostringstream ost;
			ost.str("");
			ost << m_Surfaces[SourceID]->m_shortname << "_remeshedUsingScalars";

			if (ReplaceSurface)
			{
				m_Surfaces[SourceID]->m_polyData->DeepCopy(pd);
				m_Surfaces[SourceID]->m_shortname = ost.str();
				m_Surfaces[SourceID]->UpdateScalarProperties();
			}
			else
			{
				CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
				surfProbs->m_shortname = ost.str();

				vtkTransform *sourceTrans = vtkTransform::New();
				sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

				vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
				ptransSource->SetInputData(pd);
				ptransSource->SetTransform(sourceTrans);
				ptransSource->Update();

				surfProbs->InitialiseSurface(ptransSource->GetOutput());

				AddSurfaceToRenderer(surfProbs);

				ptransSource->Delete();
				sourceTrans->Delete();
			}
		}

		pd->Delete();
		pdLengths->Delete();
	}
	else
	{
		vtkPolyData *pd = vtkPolyData::New();

		CGELRemeshing remesher;
		bool res = remesher.RemeshDirectMesh(m_Surfaces[SourceID]->m_polyData, pd, UniformTargetLength);

		if (!res)
		{
			// TODO: Update QMessageBox::information(this, "Warning","Something went wrong in remeshing");
		}
		else
		{
			std::ostringstream ost;
			ost.str("");
			ost << m_Surfaces[SourceID]->m_shortname << "_UniformRemeshed";

			if (ReplaceSurface)
			{
				m_Surfaces[SourceID]->m_polyData->DeepCopy(pd);
				m_Surfaces[SourceID]->m_shortname = ost.str();
				m_Surfaces[SourceID]->UpdateScalarProperties();
			}
			else
			{
				CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
				surfProbs->m_shortname = ost.str();

				vtkTransform *sourceTrans = vtkTransform::New();
				sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

				vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
				ptransSource->SetInputData(pd);
				ptransSource->SetTransform(sourceTrans);
				ptransSource->Update();

				surfProbs->InitialiseSurface(ptransSource->GetOutput());

				AddSurfaceToRenderer(surfProbs);

				ptransSource->Delete();
				sourceTrans->Delete();
			}
		}
		pd->Delete();
	}
}


bool C3DScene::DoICPAlignment( CICPParameters& parms )
{
	if (parms.SourceID >= m_Surfaces.size() || 
		parms.TargetID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected objects not valid!");
		return false;
	}

	if (parms.ApplyTransformTo == CICPParameters::eSource ||
		parms.ApplyTransformTo == CICPParameters::eSourceCopy)
	{
		parms.transSurfaceID = parms.SourceID;
	}

	parms.source = m_Surfaces[parms.SourceID]->m_polyData;
	parms.target = m_Surfaces[parms.TargetID]->m_polyData;
	parms.move   = m_Surfaces[parms.transSurfaceID]->m_polyData;

	if (parms.target->GetNumberOfPolys() < 1)
	{
		// TODO: Update QMessageBox::information(this, "Warning","Target must be a triangulated surface");
		return false;
	}

	parms.sourcePreTransform = m_Surfaces[parms.SourceID]->m_actor->GetMatrix();
	parms.targetPreTransform = m_Surfaces[parms.TargetID]->m_actor->GetMatrix();
	parms.moveSurfacePreTransform = m_Surfaces[parms.transSurfaceID]->m_actor->GetMatrix();

	parms.transformedSurface = vtkPolyData::New();
//	parms.surfaceErrors = vtkPolyData::New();

	// Check for normals. If signed compare and no normals for the target, then compute them
	vtkDataArray *normals = m_Surfaces[parms.TargetID]->m_polyData->GetPointData()->GetNormals(); 
	if (normals)
	{
		CProcess3DData::DoICP(parms);
	}
	else
	{
		vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
		norms->SetInputData(m_Surfaces[parms.TargetID]->m_polyData);
		norms->SplittingOff();
		norms->ConsistencyOn();
		norms->Update();

		parms.target = norms->GetOutput();

		CProcess3DData::DoICP(parms);
		norms->Delete();
	}

	std::string name;
	
//	if (parms.DoICP)
		name = m_Surfaces[parms.transSurfaceID]->m_shortname + "_ICPAlignedTo_" +
						m_Surfaces[parms.TargetID]->m_shortname;
//	else
//		name = m_Surfaces[parms.SourceID]->m_shortname + "_ComparedTo_"
//						+ m_Surfaces[parms.TargetID]->m_shortname;

	if (parms.ApplyTransformTo == CICPParameters::eSource)
	{
		m_Surfaces[parms.SourceID]->m_polyData->DeepCopy(parms.transformedSurface);
		m_Surfaces[parms.SourceID]->m_shortname = name;
		m_Surfaces[parms.SourceID]->UpdateScalarProperties();

		// Now reset source actor transform
		m_Surfaces[parms.SourceID]->m_actor->SetPosition(0,0,0);
		m_Surfaces[parms.SourceID]->m_actor->SetScale(1);
		m_Surfaces[parms.SourceID]->m_actor->SetOrientation(0,0,0);

		SetScalarBarVisible(true);
	}
	else
	{
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs2->InitialiseSurface(parms.transformedSurface);
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}
	parms.transformedSurface->Delete();
	return true;
}

void C3DScene::FindClosestPoints(unsigned int SourceID, unsigned int TargetID,bool ReplaceSurface, double MaximumDistance)
{
	if (SourceID >= m_Surfaces.size() || 
		TargetID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkMatrix4x4 *sourceMat = m_Surfaces[SourceID]->m_actor->GetMatrix();
	vtkMatrix4x4 *targetMat = m_Surfaces[TargetID]->m_actor->GetMatrix();

	vtkTransform *sourceTrans = vtkTransform::New();
	 sourceTrans->SetMatrix(sourceMat);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	 ptransSource->SetInputData(m_Surfaces[SourceID]->m_polyData);
	 ptransSource->SetTransform(sourceTrans);
	 ptransSource->Update();

	vtkTransform *targetTrans = vtkTransform::New();
	 targetTrans->SetMatrix(targetMat);

	vtkTransformPolyDataFilter *ptransTarget = vtkTransformPolyDataFilter::New();
	 ptransTarget->SetInputData(m_Surfaces[TargetID]->m_polyData);
	 ptransTarget->SetTransform(targetTrans);
	 ptransTarget->Update();


	vtkPolyData *source = ptransSource->GetOutput();
	vtkPolyData *target = ptransTarget->GetOutput();

	if(source->GetNumberOfLines() != 0 || source->GetNumberOfStrips()!=0 || target->GetNumberOfLines() != 0 || target->GetNumberOfStrips()!=0)
	{
		// TODO: Update QMessageBox::information(this, "Warning","Lines and Strips will be lost in the resulting object");
	}

	vtkPoints *npoints = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);


	vtkPointLocator *CloudLocator = vtkPointLocator::New();

	CloudLocator->SetDataSet(source);
	CloudLocator->BuildLocator();

	double pos[3];
	double psource[3];
	double actualDist;
	double val;
	vtkIdType p;
	vtkIdType id;

	int *newIndex = new int[source->GetNumberOfPoints()];
	for(int i = 0; i < source->GetNumberOfPoints(); i++)
	{
	   newIndex[i] = -1;
	}
	int nextId = 0;

    vtkDoubleArray *sca = vtkDoubleArray::SafeDownCast(source->GetPointData()->GetScalars());
	for(int j = 0; j < target->GetNumberOfPoints(); j++)
	{
		target->GetPoint(j, pos);
		p = CloudLocator->FindClosestPointWithinRadius(MaximumDistance, pos, actualDist);

		// A point was located
		if(p != -1)
		{
		  // The point was not located previously
		  if(newIndex[p] == -1)
		  {
		    newIndex[p] = nextId++;

		    source->GetPoint(p, psource);

		    verts->InsertNextCell(1);
		    id = npoints->InsertNextPoint(psource);
		    verts->InsertCellPoint(id);

			// Save scalar value if any
		    if (sca != NULL)
			{
			  val = sca->GetValue(p);
			  scalars->InsertNextTuple(&val);
			}
		  }
		}
	}

	vtkPolyData *pd = vtkPolyData::New();

	pd->SetPoints(npoints);
	npoints->Delete();
	if(sca != NULL)
		pd->GetPointData()->SetScalars(scalars);

	if(source->GetNumberOfPolys() > 0) 
	{
	  vtkCellArray *newPolys = vtkCellArray::New();
	  	  
	  vtkIdType npts;
//	  vtkIdType *pts = new vtkIdType[];
	  const vtkIdType *pts = NULL;
	  bool keep;
	  source->GetPolys()->InitTraversal();
	  while(source->GetPolys()->GetNextCell(npts, pts ))
	  {
  	     // check that all vertex are present in the resulting object
		 keep = true;
		 vtkIdList *IdL = vtkIdList::New();
		 for(int j = 0; j < npts; j++)
		 {
		    if(newIndex[pts[j]] == -1)	
			{
			  keep = false;
			  break;
			}
			else
			{
			  IdL->InsertNextId(newIndex[pts[j]]);
			}
		 }
		 if(keep)
		 { // insert polygon in the resulting object
			newPolys->InsertNextCell(IdL);
		 }
		 IdL->Delete();
	  }
	  pd->SetPolys(newPolys);
	  newPolys->Delete();
    }

	pd->SetVerts(verts);
	verts->Delete();

	if(pd->GetNumberOfPoints() == 0)
	{
	  // TODO: Update QMessageBox::information(this, "Warning","Find Closest point resulted in zero points, new point cloud was not created!");
	}
	else
	{
 	  if (ReplaceSurface)
	  {
		m_Surfaces[SourceID]->m_polyData->DeepCopy(pd);
	  }
	  else
	  {
		std::string name = m_Surfaces[SourceID]->m_shortname + "_ClosestPointsTo_" + m_Surfaces[TargetID]->m_shortname;
		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs->InitialiseSurface(pd);
		surfProbs->m_shortname = name;

		AddSurfaceToRenderer(surfProbs);
	  }
	}
	pd->Delete();
}


void C3DScene::CreateCube(double XLength, double YLength, double ZLength)
{
	vtkCubeSource *Cube = vtkCubeSource::New();
	 Cube->SetXLength(XLength);
	 Cube->SetYLength(YLength);
	 Cube->SetZLength(ZLength);
	 Cube->Update();

	std::ostringstream name;
	name << "Cube_x_" << XLength << "_y_" << YLength << "_z_" << ZLength;

	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs->InitialiseSurface(Cube->GetOutput());
	surfProbs->m_shortname = name.str();

//	m_Surfaces.push_back(surfProbs);
//	m_Renderer->AddActor(surfProbs->m_actor);
	AddSurfaceToRenderer(surfProbs);

	Cube->Delete();
}

void C3DScene::CreateTrapezoid(double XLength, double YLength, double ZLength, double alpha, double beta, bool capping)
{
	vtkTrapezoidSource *Trapezoid = vtkTrapezoidSource::New();
	 Trapezoid->SetXLength(XLength);
	 Trapezoid->SetYLength(YLength);
	 Trapezoid->SetZLength(ZLength);
	 Trapezoid->SetAlpha(alpha);
	 Trapezoid->SetBeta(beta);
	 Trapezoid->SetCapping(capping);
	 Trapezoid->Update();

	std::ostringstream name;
	name << "Trapezoid_x_" << XLength << "_y_" << YLength << "_z_" << ZLength << "_alpha_" << alpha << "_beta_" << beta << "_cap_" << capping;

	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs->InitialiseSurface(Trapezoid->GetOutput());
	surfProbs->m_shortname = name.str();

//	m_Surfaces.push_back(surfProbs);
//	m_Renderer->AddActor(surfProbs->m_actor);
	AddSurfaceToRenderer(surfProbs);

	Trapezoid->Delete();
}


void C3DScene::CreateSphere(double radius, double startTheta, double endTheta, double startPhi, double endPhi,
			int thetaRes, int phiRes, double *pos)
{
	vtkSphereSource *Sphere = vtkSphereSource::New();
	 Sphere->SetRadius(radius);
	 Sphere->SetStartTheta(startTheta);
	 Sphere->SetEndTheta(endTheta);
	 Sphere->SetStartPhi(startPhi);
	 Sphere->SetEndPhi(endPhi);
	 Sphere->SetThetaResolution(thetaRes);
	 Sphere->SetPhiResolution(phiRes);
	 Sphere->SetCenter(pos);
	 Sphere->Update();

	std::ostringstream name;
	name << "Sphere_r_" << radius << "_sT_" << startTheta << "_eT_" << endTheta
		<< "_sP_" << startPhi << "_eP_" << endPhi << "_tR_" << thetaRes << "_pR_" << phiRes;

	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs->InitialiseSurface(Sphere->GetOutput());
	surfProbs->m_shortname = name.str();

	AddSurfaceToRenderer(surfProbs);

	Sphere->Delete();
}


void C3DScene::CreateCylinder(double radius, double height, int resolution, bool capping)
{
	vtkCylinderSource *cylinder = vtkCylinderSource::New();
	 cylinder->SetRadius(radius);
	 cylinder->SetHeight(height);
	 cylinder->SetCapping(capping);
	 cylinder->SetResolution(resolution);
	 cylinder->Update();

	std::ostringstream name;
	name << "Cylinder_r_" << radius << "_h_" << height << "_res_" << resolution << "_cap_" << capping;

	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs->InitialiseSurface(cylinder->GetOutput());
	surfProbs->m_shortname = name.str();

	AddSurfaceToRenderer(surfProbs);
//	m_Surfaces.push_back(surfProbs);
//	m_Renderer->AddActor(surfProbs->m_actor);

	cylinder->Delete();
}

void C3DScene::SetupScalarBar()
{
	vtkTextProperty *SBTextProps = vtkTextProperty::New();
	 SBTextProps->ShadowOff();
	 SBTextProps->SetColor(0.1, 0.1, 0.1);
//	 SBTextProps->BoldOn();

	m_scalarBar = vtkScalarBarActor::New();
	 m_scalarBar->SetLookupTable(m_lookup);
	 m_scalarBar->SetPosition(0.01 , 0.1 );
//	 scalarBar->SetPosition2(0.42 , 0.7 );
//	 scalarBar->SetOrientationToHorizontal();
	 m_scalarBar->SetOrientationToVertical();
	 m_scalarBar->SetWidth(m_scalarBar->GetWidth() / 2);
	 m_scalarBar->SetLabelTextProperty(SBTextProps);
	 m_scalarBar->SetVisibility(m_ScalarBarVisible);
	 m_scalarBar->SetMaximumNumberOfColors(255);

	SBTextProps->Delete();

	m_Renderer->AddActor(m_scalarBar);
}



void C3DScene::SetScalarBarVisible(bool flag)
{
	if (m_ScalarBarVisible != flag)
	{
		m_ScalarBarVisible = flag;

		m_scalarBar->SetVisibility(m_ScalarBarVisible);
	}
}

bool C3DScene::GetScalarBarVisible() const
{
	return m_ScalarBarVisible;
}

void C3DScene::SetAxesVisible(bool flag)
{
	if (m_AxesVisible != flag)
	{
		m_AxesVisible = flag;

		m_Axes->SetVisibility(m_AxesVisible);

		if (m_AxesVisible)
		{			
			m_Renderer->AddActor(m_Axes);
		}
		else
		{
			m_Renderer->RemoveActor(m_Axes);
		}
	}
}

bool C3DScene::GetAxesVisible() const
{
	return m_AxesVisible;
}



std::string C3DScene::GetSurfaceShortName(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return "";

	return m_Surfaces[id]->m_shortname;
}

std::string C3DScene::GetSurfaceFullName(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return "";

	return m_Surfaces[id]->m_fullname;
}

bool C3DScene::GetUndoAvaible(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return false;

	return m_Surfaces[id]->UndoAvailable();
}


vtkActor * C3DScene::GetActor(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return NULL;

	return m_Surfaces[id]->m_actor;
}

vtkMapper * C3DScene::GetMapper(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return NULL;

	return m_Surfaces[id]->m_mapper;
}

bool C3DScene::GetSurfaceScalarRange(unsigned int id, double *range)
{
	if (id >= m_Surfaces.size())
		return false;

	if (m_Surfaces[id]->m_polyData->GetPointData()->GetScalars() == NULL)
		return false;

	m_Surfaces[id]->m_polyData->GetPointData()->GetScalars()->GetRange(range);
	return true;
}

bool C3DScene::SurfaceHasNormals(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return false;

	vtkDataArray *normals = m_Surfaces[id]->m_polyData->GetPointData()->GetNormals(); 

	if (normals)
		return true;

	return false;
}

std::string C3DScene::GetSurfaceValues(unsigned int id) const
{
	if (id >= m_Surfaces.size())
		return "Surface ID out of range";

	vtkPolyData *pd = m_Surfaces[id]->m_polyData;
	return vtkExtMisc::GetSurfaceValues(pd, "\r\n");
}

void C3DScene::SetScalarRange(unsigned int id, double *range)
{
	if (id >= m_Surfaces.size())
		return;

	m_Surfaces[id]->SetScalarRange(range);
}

bool C3DScene::GetScalarRange(unsigned int id, double *range)
{
	if (id >= m_Surfaces.size())
		return false;

	return m_Surfaces[id]->GetScalarRange(range);
}

bool C3DScene::SavePlane(const std::string& fname) const
{
	if (m_PlaneWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No plane!");
		return false;
	}

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkPolyData *pd = vtkPolyData::New();

	vtkIdType id = pts->InsertNextPoint(m_PlaneWidget->GetOrigin());
	verts->InsertNextCell(1);
	verts->InsertCellPoint(id);

	id = pts->InsertNextPoint(m_PlaneWidget->GetPoint1());
	verts->InsertNextCell(1);
	verts->InsertCellPoint(id);

	id = pts->InsertNextPoint(m_PlaneWidget->GetPoint2());
	verts->InsertNextCell(1);
	verts->InsertCellPoint(id);

	pd->SetPoints(pts);
	pd->SetVerts(verts);

	bool result = vtkExtMisc::WritePDVTK(pd, fname);
	
	pts->Delete();
	verts->Delete();
	pd->Delete();

	return result;
}

bool C3DScene::LoadPlane(const std::string& fname)
{
	if (m_PlaneWidget == NULL)
	{
		// TODO: Update QMessageBox::information(this, "Warning","No plane!");
		return false;
	}

	vtkPolyDataReader *pd = vtkExtMisc::SafeReadPolyData(fname);
	if (!pd) return false;

	const int N = pd->GetOutput()->GetNumberOfPoints();

	if (N != 3)
		return false;

	m_PlaneWidget->SetOrigin(pd->GetOutput()->GetPoint(0));
	m_PlaneWidget->SetPoint1(pd->GetOutput()->GetPoint(1));
	m_PlaneWidget->SetPoint2(pd->GetOutput()->GetPoint(2));

	pd->Delete();
	return true;
}

void C3DScene::UpdateAllBounds()
{
	if (m_Surfaces.size() > 0)
	{
		for (int j = 0; j < 6; j++)
				m_AllBounds[j] = 0;
		
		for (unsigned int i = 0; i < m_Surfaces.size(); i++)
		{
			double b[6];
			m_Surfaces[i]->m_polyData->GetBounds(b);

			m_AllBounds[0] = std::min(m_AllBounds[0], b[0]);
			m_AllBounds[1] = std::max(m_AllBounds[1], b[1]);
			m_AllBounds[2] = std::min(m_AllBounds[2], b[2]);
			m_AllBounds[3] = std::max(m_AllBounds[3], b[3]);
			m_AllBounds[4] = std::min(m_AllBounds[4], b[4]);
			m_AllBounds[5] = std::max(m_AllBounds[5], b[5]);
		}
		double xl = m_AllBounds[1] - m_AllBounds[0];
		double yl = m_AllBounds[3] - m_AllBounds[2];
		double zl = m_AllBounds[5] - m_AllBounds[4];

		double l = std::max(std::max(xl, yl), zl);

		m_Axes->SetTotalLength(l / 2, l / 2, l / 2);
	}
}




void C3DScene::SetMarkerValue(double val)
{
	m_MarkerValue = val;
}

double C3DScene::GetMarkerValue() const
{
	return m_MarkerValue;
}

void C3DScene::SetStatusText(const std::string txt, bool visible)
{
	m_StatusText->SetInput(txt.c_str());
	m_StatusText->SetVisibility(visible);
}

void C3DScene::SetGradientBackground(bool flag)
{
	if (m_BackgroundActor)
	{
		m_BackgroundActor->SetVisibility(flag);
	}
}

bool C3DScene::GetGradientBackground() const
{
	if (m_BackgroundActor)
	{
		return m_BackgroundActor->GetVisibility() == 1;
	}
	return false;
}

void C3DScene::RotateCamera()
{
	m_Renderer->GetActiveCamera()->Azimuth(1);
}

void C3DScene::Animate( unsigned int SourceID, unsigned int DisplacementID, unsigned int &TargetID, double t )
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected objects not valid!");
		return;
	}
	if (TargetID > m_Surfaces.size())
	{
		// Create new output surface
		std::string name = m_Surfaces[SourceID]->m_shortname + "_animated";
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);

		surfProbs2->InitialiseSurface(m_Surfaces[SourceID]->m_polyData);
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);

		TargetID = m_Surfaces.size() - 1;
		return;
	}

	if (TargetID >= m_Surfaces.size() ||  DisplacementID >= m_Surfaces.size())
		return;


	int N = m_Surfaces[SourceID]->m_polyData->GetNumberOfPoints();
	if (N != m_Surfaces[DisplacementID]->m_polyData->GetNumberOfPoints() || N != m_Surfaces[TargetID]->m_polyData->GetNumberOfPoints())
		return;

	for (int i = 0; i < N;  i++)
	{
		double p[3];
		double d[3];

		m_Surfaces[SourceID]->m_polyData->GetPoint(i, p);
		m_Surfaces[DisplacementID]->m_polyData->GetPoint(i, d);

		p[0] = p[0] + t * d[0];
		p[1] = p[1] + t * d[1];
		p[2] = p[2] + t * d[2];

		m_Surfaces[TargetID]->m_polyData->GetPoints()->SetPoint(i, p);
	}
	m_Surfaces[TargetID]->m_polyData->Modified();
}

void C3DScene::DoProjectSurface( unsigned int SourceID, unsigned int TargetID, bool CreateDisplacement, bool ReplaceSource )
{
	if (SourceID >= m_Surfaces.size() || 
		TargetID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning", "Selected objects not valid!");
		return;
	}

	vtkPolyData *source = m_Surfaces[SourceID]->m_polyData;
	vtkPolyData *target = m_Surfaces[TargetID]->m_polyData;
	vtkMatrix4x4 *sourcePreTransform = m_Surfaces[SourceID]->m_actor->GetMatrix();
	vtkMatrix4x4 *targetPreTransform = m_Surfaces[TargetID]->m_actor->GetMatrix();

	// Check for normals. If signed compare and no normals for the target, then compute them
	vtkDataArray *normals = target->GetPointData()->GetNormals(); 

	vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
	if (!normals)
	{
		norms->SetInputData(target);
		norms->SplittingOff();
		norms->ConsistencyOn();
		norms->Update();

		target = norms->GetOutput();
	}

	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(sourcePreTransform);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(source);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	vtkTransform *targetTrans = vtkTransform::New();
	targetTrans->SetMatrix(targetPreTransform);

	vtkTransformPolyDataFilter *ptransTarget = vtkTransformPolyDataFilter::New();
	ptransTarget->SetInputData(target);
	ptransTarget->SetTransform(targetTrans);
	ptransTarget->Update();

	vtkPolyDataProjection *PDProj = vtkPolyDataProjection::New();
	PDProj->SetCreateDisplacementField(CreateDisplacement);

	PDProj->SetInputConnection(ptransSource->GetOutputPort());
	PDProj->SetTargetData(ptransTarget->GetOutput());
	PDProj->Update();

	std::string name = m_Surfaces[SourceID]->m_shortname + "_ProjectedTo_"	+ m_Surfaces[TargetID]->m_shortname;

	if (ReplaceSource)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(PDProj->GetOutput());
		m_Surfaces[SourceID]->m_shortname = name;
		m_Surfaces[SourceID]->UpdateScalarProperties();
		SetScalarBarVisible(true);
	}
	else
	{
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs2->InitialiseSurface(PDProj->GetOutput());
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}
	if (CreateDisplacement)
	{
		name = m_Surfaces[SourceID]->m_shortname + "_ProjectedTo_"	+ m_Surfaces[TargetID]->m_shortname + "_displacements";

		CSurfaceProperties *surfProbs3 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs3->InitialiseSurface(PDProj->GetDisplacementField());
		surfProbs3->m_shortname = name;

		AddSurfaceToRenderer(surfProbs3);
	}

	PDProj->Delete();
	ptransTarget->Delete();
	targetTrans->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
	norms->Delete();
}

void C3DScene::RemoveNormals( unsigned int id )
{
	if (id >= m_Surfaces.size())
		return;

	m_Surfaces[id]->m_polyData->GetPointData()->SetNormals(NULL);
}

void C3DScene::RemoveScalars( unsigned int id )
{
	if (id >= m_Surfaces.size())
		return;

	m_Surfaces[id]->m_polyData->GetPointData()->SetScalars(NULL);
	
	double range[2];
	range[0] = 0; range[0] = 0;
	m_Surfaces[id]->SetScalarRange(range);
}

void C3DScene::CalculatePCANormalsParameters( unsigned int SourceID, double RadiusFactor,  
											 double &ManualNormalRadius, double &ManualPlaneDist, 
											 int &ManualMinCliqueSize, double &CliqueNeigDist )
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	//CMRFSurfaceReconstruction::CParms parms;
	//parms.SetInputType(1);
	//parms.EstimateNormals = true;
	//parms.WriteLevel = 0;

	//parms.AdaptiveParameters = true;
	//parms.NormalRadiusFactor = RadiusFactor;

	//std::string msg;
	//CMRFSurfaceReconstruction MRFSurface;
	//bool result = MRFSurface.AdaptiveParameters(m_Surfaces[SourceID]->m_polyData, parms);
	//if (!result)
	//{
	//	// TODO: Update QMessageBox::information(this, "Warning","Could not compute adaptive parameters");
	//}
	//ManualNormalRadius = parms.NormalRadius;
	//ManualPlaneDist = parms.MaxPlaneDistance;
	//ManualMinCliqueSize = parms.MinCliqueSize;
	//CliqueNeigDist = parms.CliqueNeighbourDistance;
}

bool C3DScene::ComputeAverageEdgelengths(unsigned int SourceID, double & avgL)
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return false;
	}

	avgL = CGELRemeshing::ComputeAverageEdgeLength(m_Surfaces[SourceID]->m_polyData);
	if (avgL == 0)
		return false;

	return true;
}

void C3DScene::CalculatePCANormals( unsigned int SourceID, bool ReplaceSurface, bool ConComp, bool KeepLargestClique,
								   bool UseOrgNormals, bool VisNormals, int AdaptiveOrManual, double RadiusFactor, 
								   double PointsPerNormal, double ManualNormalRadius, double ManualPlaneDist, 
								   int ManualMinCliqueSize, double CliqueNeigDist )
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	//vtkTransform *sourceTrans = vtkTransform::New();
	//sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

	//vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	//ptransSource->SetInputData(m_Surfaces[SourceID]->m_polyData);
	//ptransSource->SetTransform(sourceTrans);
	//ptransSource->Update();


	//CMRFSurfaceReconstruction::CParms parms;
	//parms.SetInputType(1);
	//parms.EstimateNormals = true;
	//parms.ConnectedComponentAnalysis = ConComp;
	//parms.LargestCliqueOnly = KeepLargestClique;
	//parms.UseReferenceNormalsInVoting = UseOrgNormals;
	//parms.WriteLevel = 0;


	//if (AdaptiveOrManual == 0)
	//{
	//	// Adaptive parameters
	//	parms.AdaptiveParameters = true;
	//	parms.NormalRadiusFactor = RadiusFactor;
	//	parms.NormalSearchMode = 0;
	//	if (PointsPerNormal > 0)
	//	{
	//		parms.NormalSearchMode = 1;
	//		parms.NormalPointsPerNormal	 = PointsPerNormal;
	//	}
	//}
	//else
	//{
	//	// Manual parameters
	//	parms.AdaptiveParameters = false;
	//	parms.NormalRadius = ManualNormalRadius;
	//	parms.NormalSearchMode = 0;
	//	if (PointsPerNormal > 0)
	//	{
	//		parms.NormalSearchMode = 1;
	//		parms.NormalPointsPerNormal	 = PointsPerNormal;
	//	}
	//	parms.MinCliqueSize = ManualMinCliqueSize;
	//	parms.MaxPlaneDistance = ManualPlaneDist;
	//	parms.CliqueNeighbourDistance = CliqueNeigDist;
	//}

	//std::string msg;
	//CMRFSurfaceReconstruction MRFSurface;
	//bool result = MRFSurface.OneShotComputePCANormals(parms,ptransSource->GetOutput(), msg);
	//if (!result)
	//{
	//	// TODO: Update QMessageBox::information(this, "Warning",msg.c_str());
	//}
	//else
	//{
	//	if (ReplaceSurface)
	//	{
	//		m_Surfaces[SourceID]->m_polyData->DeepCopy(MRFSurface.GetNormalData());
	//	}
	//	else
	//	{
	//		std::string name = m_Surfaces[SourceID]->m_shortname + "_PCANormals";
	//		CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	//		surfProbs->InitialiseSurface(MRFSurface.GetNormalData());
	//		surfProbs->m_shortname = name;
	//		AddSurfaceToRenderer(surfProbs);
	//	}
	//	if (VisNormals)
	//	{
	//		if (ReplaceSurface)
	//		{
	//			VisualiseNormals(SourceID, false);
	//		}
	//		else
	//		{
	//			VisualiseNormals(m_Surfaces.size()-1, false);
	//		}
	//	}
	//}

	//ptransSource->Delete();
	//sourceTrans->Delete();
}

std::string C3DScene::PickPointAndGetText( double x, double y )
{
	std::ostringstream ost;
	std::string pointInf;

	int res = m_PointPicker->Pick(x, y, 0, m_Renderer);

	bool found = false;
	if (res)
	{
		int id = m_PointPicker->GetPointId();

		vtkActor *pact = m_PointPicker->GetActor();

		if (pact != NULL && id >= 0)
		{
			double p[3];
			pact->GetMapper()->GetInput()->GetPoint(id, p);

			m_PickSphereActor->SetPosition(p);
			m_PickSphereActor->SetVisibility(1);

			ost << "Point ID " << id << " (" << p[0] << ", " << p[1] << ", " << p[2] << ")";

			vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pact->GetMapper()->GetInput()->GetPointData()->GetScalars());
			if (!scalars)
			{
				if (pact->GetMapper()->GetInput()->GetPointData()->GetScalars() && pact->GetMapper()->GetInput()->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
				{
					vtkUnsignedCharArray *scalU  = vtkUnsignedCharArray::SafeDownCast(pact->GetMapper()->GetInput()->GetPointData()->GetScalars());
					if (scalU)
					{
						unsigned char val[3];
						scalU->GetTypedTuple(id, val);
						ost << " value (" << (int)val[0] << ", " << (int)val[1] << ", " << (int)val[2] << ")";
					}
					else
					{
						ost << " 3-component scalar ";
					}
				}
				else if (pact->GetMapper()->GetInput()->GetPointData()->GetScalars())
				{
					ost << pact->GetMapper()->GetInput()->GetPointData()->GetScalars()->GetNumberOfComponents() << " component scalars";
				}
			}
			else
			{
				double scalar = pact->GetMapper()->GetInput()->GetPointData()->GetScalars()->GetTuple1(id);
				ost << " value " << scalar;
			}


			if (pact->GetMapper()->GetInput()->GetPointData()->GetTCoords())
			{
				double *tc = pact->GetMapper()->GetInput()->GetPointData()->GetTCoords()->GetTuple(id);
				ost << " TCoords (" << tc[0] << ", " << tc[1] << ")";
			}

//			ost << "Points ID " << id << " (" << p[0] << ", " << p[1] << ", " << p[2] << ") POINT (" << x << ", " << y << ")";

			found = true;

			// Resize sphere
			double bounds[6];
			pact->GetMapper()->GetInput()->GetBounds(bounds);

			double l = sqrt((bounds[1]-bounds[0]) * (bounds[1]-bounds[0]) + 
				(bounds[3]-bounds[2]) * (bounds[3]-bounds[2]) +
				(bounds[5]-bounds[4]) * (bounds[5]-bounds[4]));

			// Sphere radius 0.1 % of bounding box diagonal
			double NL = l * 0.001;

			m_PickSphere->SetRadius(NL);
		}
	}
	
	if (!found)
	{
		//m_PickSphereActor->SetVisibility(0);
		//ost << "No pick found when POINT (" << x << ", " << y << ")";
	}

	pointInf = ost.str();

	return pointInf;
}

void C3DScene::SubSamplePointCloud( unsigned int SourceID, bool ReplaceSurface, int NumPoints )
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkPolyData *source = m_Surfaces[SourceID]->m_polyData;
	vtkMatrix4x4 *sourcePreTransform = m_Surfaces[SourceID]->m_actor->GetMatrix();

	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(sourcePreTransform);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(source);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	vtkPolyDataResamplePoints *Resample = vtkPolyDataResamplePoints::New();
	Resample->SetInputConnection(ptransSource->GetOutputPort());
	Resample->SetNumberOfOutputPoints(NumPoints);
	Resample->Update();

	std::string name = m_Surfaces[SourceID]->m_shortname + "_PointSubSampled";

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(Resample->GetOutput());
		m_Surfaces[SourceID]->m_shortname = name;
		m_Surfaces[SourceID]->UpdateScalarProperties();
		SetScalarBarVisible(true);
	}
	else
	{
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs2->InitialiseSurface(Resample->GetOutput());
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}

	Resample->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();

}


void C3DScene::RandomResample( unsigned int SourceID, bool ReplaceSurface, int NumPoints, bool CloneNormals )
{

	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}

	vtkPolyData *source = m_Surfaces[SourceID]->m_polyData;
	vtkMatrix4x4 *sourcePreTransform = m_Surfaces[SourceID]->m_actor->GetMatrix();

	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(sourcePreTransform);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(source);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	vtkPolyDataRandomResamplePoint *randomResample = vtkPolyDataRandomResamplePoint::New();
	randomResample->SetInputConnection(ptransSource->GetOutputPort());
	randomResample->SetResamplePoints(NumPoints);
	randomResample->SetCreatePointNormals(CloneNormals);
	randomResample->Update();

	std::string name = m_Surfaces[SourceID]->m_shortname + "_PointResampled";

	if (ReplaceSurface)
	{
		m_Surfaces[SourceID]->m_polyData->DeepCopy(randomResample->GetOutput());
		m_Surfaces[SourceID]->m_shortname = name;
		m_Surfaces[SourceID]->UpdateScalarProperties();
		SetScalarBarVisible(true);
	}
	else
	{
		CSurfaceProperties *surfProbs2 = new CSurfaceProperties(m_lookup, mSettings);
		surfProbs2->InitialiseSurface(randomResample->GetOutput());
		surfProbs2->m_shortname = name;

		AddSurfaceToRenderer(surfProbs2);
	}

	randomResample->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
}

void C3DScene::NormalsFromCamera( unsigned int SourceID )
{
	if (SourceID >= m_Surfaces.size())
	{
		// TODO: Update QMessageBox::information(this, "Warning","Selected object not valid!");
		return;
	}
	
	double CameraPos[3];
	m_Renderer->GetActiveCamera()->GetPosition(CameraPos);

	vtkDoubleArray *Normals = vtkDoubleArray::New();
	Normals->SetNumberOfComponents(3);
	Normals->SetName("Normals");

	// Use actor transform before calculating normals
	vtkTransform *sourceTrans = vtkTransform::New();
	sourceTrans->SetMatrix(m_Surfaces[SourceID]->m_actor->GetMatrix());

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	ptransSource->SetInputData(m_Surfaces[SourceID]->m_polyData);
	ptransSource->SetTransform(sourceTrans);
	ptransSource->Update();

	for (int i = 0; i < ptransSource->GetOutput()->GetNumberOfPoints(); i++)
	{
		double p[3];
		double n[3];

		ptransSource->GetOutput()->GetPoint(i, p);
		n[0] = CameraPos[0] - p[0];
		n[1] = CameraPos[1] - p[1];
		n[2] = CameraPos[2] - p[2];

		vtkMath::Normalize(n);
		Normals->InsertNextTuple(n);
	}
	m_Surfaces[SourceID]->m_polyData->GetPointData()->SetNormals(Normals);
	Normals->Delete();

	ptransSource->Delete();
	sourceTrans->Delete();
}

void C3DScene::SetScalarLookupTableNum( int num )
{
	//m_lookup->SetNumberOfColors(5);
	//m_lookup->SetTableValue(0, 0, 0, 1, 1);
	//m_lookup->SetTableValue(1, 0, 1, 0, 1);
	//m_lookup->SetTableValue(2, 0, 1, 1, 1);
	//m_lookup->SetTableValue(3, 1, 0, 0, 1);
	//m_lookup->SetTableValue(4, 1, 0, 1, 1);

	//	m_lookup->SetHueRange(0.2, 0);


	// Grey to skin
	//m_lookup->SetHueRange(40/360, 40/360);
	//m_lookup->SetSaturationRange(0, 0.30);
	//m_lookup->SetValueRange(0.93,0.93);

	// Red to skin
	//m_lookup->SetHueRange(0/360, 40/360);
	//m_lookup->SetSaturationRange(1.0, 0.30);
	//m_lookup->SetValueRange(1.0,0.93);

	// blue to skin
	//m_lookup->SetHueRange(240.0/360.0, 40/360);
	//m_lookup->SetSaturationRange(1.0, 0.30);
	//m_lookup->SetValueRange(1.0,0.93);


	if (num == 0)
	{
		// Standard blue to red (rainbow)
		m_lookup->SetHueRange(2.0/3.0, 0);
		m_lookup->SetSaturationRange(1, 1);
		m_lookup->SetValueRange(1.0,1);
	}
	else if (num == 1)
	{
		// blue to skin
		m_lookup->SetHueRange(240.0/360.0, 40/360);
		m_lookup->SetSaturationRange(0.7, 0.30);
		m_lookup->SetValueRange(1.0,0.93);
	}
	else if (num == 2)
	{
		// Grey to blue
		m_lookup->SetHueRange(0.3, 0.6);
		m_lookup->SetSaturationRange(0, 0.8);
		m_lookup->SetValueRange(0.5,0.8);
	}
	if (num == 3)
	{
		// Inverse rainbow
		m_lookup->SetHueRange(0, 2.0/3.0);
		m_lookup->SetSaturationRange(1, 1);
		m_lookup->SetValueRange(1.0,1);
	}
	if (num == 4)
	{
		// Dark Grey to golden
		m_lookup->SetHueRange(50.0/360.0, 50.0/360.0);
		m_lookup->SetSaturationRange(0, 0.85);
		m_lookup->SetValueRange(0.5,0.8);
	}
	if (num == 5)
	{
		// Hot
		m_lookup->SetHueRange(0.0/360.0, 60.0/360.0);
		m_lookup->SetSaturationRange(1, 1);
		m_lookup->SetValueRange(0.5,1);
	}
	if (num == 6)
	{
		// Hot
		m_lookup->SetHueRange(0.0/360.0, 0.0/360.0);
		m_lookup->SetSaturationRange(0, 0);
		m_lookup->SetValueRange(0.2,0.8);
	}
	// Grey to bluish-purplish
	//m_lookup->SetHueRange(280.0/360.0, 280.0/360.0);
	//m_lookup->SetSaturationRange(0, 0.8);
	//m_lookup->SetValueRange(0.8,0.8);

	// Dark Grey to bluish-purplish
	//m_lookup->SetHueRange(280.0/360.0, 280.0/360.0);
	//m_lookup->SetSaturationRange(0, 0.8);
	//m_lookup->SetValueRange(0.5,0.8);




	// Purplish to golden
	//m_lookup->SetHueRange(280.0/360.0, 50.0/360.0);
	//m_lookup->SetSaturationRange(0.8, 0.85);
	//m_lookup->SetValueRange(0.5,0.8);

	m_lookup->Build();
	//m_lookup->ForceBuild();
	//m_lookup->Modified();
	//if (m_scalarBar)
	//{
	//	m_scalarBar->SetLookupTable(m_lookup);
	//	m_scalarBar->Modified();
	//}
	//for (unsigned int i = 0; i < m_Surfaces.size(); i++)
	//{
	//	m_Surfaces[i]->UpdateScalarProperties();
	//}
}
//
//void C3DScene::SetDefaultPointSize( int ps )
//{
//	m_DefaultPointSize = ps;
//}

void C3DScene::SetBackgroundColor( int R, int G, int B )
{
	m_Renderer->SetBackground(R, G, B);
}

void C3DScene::SetSphereWidgetRadius( double val )
{
	if (m_SphereWidget != NULL)
	{
		m_SphereWidget->SetRadius(val);
	}
}

double C3DScene::GetSphereWidgetRadius() const
{
	if (m_SphereWidget != NULL)
	{
		return m_SphereWidget->GetRadius();
	}
	return -1;
}

void C3DScene::MergeAllSets()
{
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkPolyData *pd = vtkPolyData::New();
	vtkDoubleArray *scalD = vtkDoubleArray::New();
	vtkUnsignedCharArray *scalU = vtkUnsignedCharArray::New();  // We need different types to handle both single scalars and RGB

	bool isScalars = true;

	int NumComp = 0;
	if (m_Surfaces[0]->m_polyData->GetPointData()->GetScalars())
	{
		NumComp = m_Surfaces[0]->m_polyData->GetPointData()->GetScalars()->GetNumberOfComponents();
		scalD->SetNumberOfComponents(NumComp);
		scalU->SetNumberOfComponents(NumComp);
	}
	else
	{
		isScalars = false;
	}

	vtkDoubleArray *Normals = vtkDoubleArray::New();
	Normals->SetNumberOfComponents(3);
	Normals->SetName("Normals");

	bool IsNormals = true;
	for (unsigned int j = 0; j < m_Surfaces.size(); j++)
	{
		vtkPolyData *surf = m_Surfaces[j]->m_polyData;

		vtkDoubleArray *norms = vtkDoubleArray::SafeDownCast(surf->GetPointData()->GetNormals()); 
		if (!norms)
		{
			IsNormals = false;
		}
		vtkDoubleArray *scalsD = NULL;
		vtkUnsignedCharArray *scalsU = NULL;
		if (NumComp == 1)
			scalsD = vtkDoubleArray::SafeDownCast(surf->GetPointData()->GetScalars());
		if (NumComp == 3)
			scalsU = vtkUnsignedCharArray::SafeDownCast(surf->GetPointData()->GetScalars());

		if (!scalsD && ! scalsU)
		{
			isScalars = false;
		}

		for (int i = 0; i < surf->GetNumberOfPoints(); i++)
		{
			double p[3];
			surf->GetPoint(i, p);

			vtkIdType id = pts->InsertNextPoint(p);
			verts->InsertNextCell(1);
			verts->InsertCellPoint(id);

			if (IsNormals)
			{
				double n[3];
				norms->GetTuple(i, n);
				Normals->InsertNextTuple(n);
			}

			if (isScalars && NumComp == 1)
			{
				double val = scalsD->GetValue(i);
				scalD->InsertNextTuple(&val);
			}
			if (isScalars && NumComp == 3)
			{
				unsigned char val[3];
				scalsU->GetTypedTuple(i, val);
				scalU->InsertNextTypedTuple(val);
			}
		}
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	if (IsNormals)
		pd->GetPointData()->SetNormals(Normals);	
	if (isScalars && NumComp == 1)
		pd->GetPointData()->SetScalars(scalD);
	if (isScalars && NumComp == 3)
		pd->GetPointData()->SetScalars(scalU);
	scalD->Delete();
	scalU->Delete();

	pts->Delete();
	verts->Delete();
	Normals->Delete();

	CSurfaceProperties *surfProbs = new CSurfaceProperties(m_lookup, mSettings);
	surfProbs->InitialiseSurface(pd);
	surfProbs->m_shortname = "AllMerged";

	AddSurfaceToRenderer(surfProbs);

	pd->Delete();
}

