#include "MRFSurfaceReconstruction.h"

#include <sys/stat.h>

#include <sstream>
#include <string>
#include <vtkMath.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkImageResample.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkImplicitModeller.h>
#include <vtkTimerLog.h>
#include <vtkImplicitVolume.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkGlyph3D.h>
#include <vtkLineSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkClipPolyData.h>
#include <vtkSTLReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkCubeSource.h>
#include <vtkTransform.h>
#include <vtkPointLocator.h>
#include <vtkSphereSource.h>
#include <vtkClipPolyData.h>
#include <vtkSphere.h>
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkImageReslice.h>
#include <vtkMetaImageWriter.h>

#include <vtkVRMLImporter.h>
#include <vtkRenderer.h>
#include <vtkMapper.h>
#include <vtkPointLocator.h>
#include <vtkFloatArray.h>
#include <vtkFeatureEdges.h>

#include "MemoryUtils.h"
#include "vtkPolyDataDifference.h"
#include "vtkOrientedPointSetDistanceFilter.h"
#include "vtkOrientedPointSetDistanceFilter2.h"
#include "vtkSignedDistanceTransformFilter.h"
#include "vtkSignedDistanceFieldLaplacianFilter.h"
#include "vtkDistanceFieldMRFRegularisation.h"
#include "BloomenthalPolygonizer.h"
#include "vtkRemoveNonManifoldCells.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkRemoveUnusedPolyDataPoints.h"
//#include "vtk3DMDTxtReader.h"
#include "vtkPolyDataPCANormals.h"
#include "vtkPolyDataOrientNormalsByVoting.h"
#include "vtkPolyDataSamplingDensity.h"
#include "vtkPrincipalAxisTransform.h"
#include "vtkPolyDataMarkCoveredCells.h"
#include "vtkCropPolyData.h"
#include "vtkPolyDataMarkUsingVolume.h"
#include "GeneralUtils.h"
#include "vtkPolyDataDifference.h"
#include "vtkExtMisc.h"
#include "GELRemeshing.h"
#include <vtkRemoveNonManifoldCells.h>
//#include <vil/vil_rgb.h>
//#include <vil/vil_load.h>
//#include <vil/vil_save.h>
//#include <vil/vil_image_view.h>



//#include <3dMDCamera.h>
#include "MultiTimer.h"


CMRFSurfaceReconstruction::CMRFSurfaceReconstruction()
{
	RawData = NULL;
	NormalData = NULL;
	MRFDistanceField = NULL;
	MRFSurface = NULL;
	MRFSurfaceCropped = NULL;
	MRFSurfaceAggresivelyCropped = NULL;
}

CMRFSurfaceReconstruction::~CMRFSurfaceReconstruction()
{
	DeleteAll();
}

void CMRFSurfaceReconstruction::DeleteAll()
{
	if (RawData)
	{
		RawData->Delete();
		RawData = NULL;
	}
	if (NormalData)
	{
		NormalData->Delete();
		NormalData = NULL;
	}
	if (MRFDistanceField)
	{
		MRFDistanceField->Delete();
		MRFDistanceField = NULL;
	}
	if (MRFSurface)
	{
		MRFSurface->Delete();
		MRFSurface = NULL;
	}
	if (MRFSurfaceCropped)
	{
		MRFSurfaceCropped->Delete();
		MRFSurfaceCropped = NULL;
	}
	if (MRFSurfaceAggresivelyCropped)
	{
		MRFSurfaceAggresivelyCropped->Delete();
		MRFSurfaceAggresivelyCropped = NULL;
	}
}

//#include <thread>
//
//void ThreadTestHello()
//{
//	std::cout << "Hello from thread" << std::endl;
//}
//
//void ThreadTest()
//{
//	std::cout << "Thread test start -----------------------------------------" << std::endl;
//	std::thread t(ThreadTestHello); 
//	t.join();
//	std::cout << "Thread test end -----------------------------------------" << std::endl;
//}
//

void CMRFSurfaceReconstruction::Compute(CMRFSurfaceReconstruction::CParms &parms)
{
//	ThreadTest();

	m_Parms = parms;
	m_Parms.InputNameShort = CGeneralUtils::StripPathAndExtensionFromFilename(m_Parms.InputName);

	if (m_Parms.WriteLevel > 1)
	{
		m_Parms.DumpToFile();
	}

//	ThreadTest();

	ReadCalibrationFiles();

	if (!ConvertRawScanData())
		return;

//	ThreadTest();

	if (!InputDataStatistics())
		return;

//	ThreadTest();

	// This is kind of hacky
	if (m_Parms.SignedDistance == false)
	{
		m_Parms.EstimateNormals = false;
		m_Parms.InputType = 5;
		m_Parms.ConnectedComponentAnalysis = false;
		m_Parms.DistanceMode = VTK_UNSIGNEDDISTANCE_MEDIAN;
		m_Parms.Remesh = 2;
	}

	if (!CreateNormals())
		return;
	
//	ThreadTest();

	if (m_Parms.UseLeavePatchOut)
	{
		CreateCutSpheres();
		CutWithSpheres();
	}
	else if (m_Parms.AddNoise)
	{
		AddNoiseToCloud();
	}

	if (m_Parms.Optimisation == VTK_OPTIMISATION_SPARSECHOL)
	{
		SparseCholeskyMRFRegularisation();
		MultiLevelMRFRefinement();
	}
	else
	{
		if (m_Parms.MultiLevel)
		{
			if (!MultiLevelMRFRegularisation())
				return;
		}
		else
		{
			CreateDensity();
			ComputeDistanceMap();
			PolygoniseRawDistanceMap();
			MRFRegularisation();
			PolygoniseDistanceMap();
			PolygoniseUpscaledDistanceMap();
		}
	}
	if (parms.DoComputeSurface)
	{
		if (!PolygoniseDistanceMap())
			return;

		//if (!CheckSurfaceConsistency())
		//	return;

		if (!CropSurface())
			return;

		if (m_Parms.WriteLevel > 1)
		{
			SurfaceStatistics();
		}


		//if (m_Parms.Remesh)
		//{
		//	if (!Remesh())
		//		return;
		//}

		//	WeightColorSurface();
		if (m_Parms.WriteLevel > 1)
		{
			ApproximateWeightColorSurface(0);
			if (m_Parms.Remesh)
			{
				ApproximateWeightColorSurface(1);
			}
		}
		if (m_Parms.WriteLevel > 1)
		{
			ComputeAccuracy();
		}

		if (m_Parms.UseLeavePatchOut)
		{
			ComputePatchAccuracy();
		}

		if (m_Parms.CloneTextures)
		{
			//		CloneTexturesFromVRML();
			ReTextureSurface();
		}
	}
}

bool CMRFSurfaceReconstruction::CreateNormals()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_raw.vtk";
	std::string outputnameNorms		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_normals.vtk";
	std::string outputnameONorms    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";

	if (NormalData)
		return true;

	if (m_Parms.InputName != "" && file_exists(outputnameONorms.c_str()))
	{
		return true;
	}

	if (!RawData)
	{
		RawData = vtkExtMisc::MultiReadSurface(iname);
		if (!RawData)
		{
			std::cerr << "Could not read " << iname << std::endl;
			return false;
		}
	}

	std::cout << "Estimating normals with input data type " << m_Parms.InputType << std::endl;	
	CMultiTimer::MultiTimer().Start("Estimating normals", 1);

	if (m_Parms.InputType == 0)
	{
		// treatments m1 and m2
		vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
		cleaner->SetInputData(RawData);
		cleaner->Update();

		int Ntot = RawData->GetNumberOfPoints();
		int NCleaned = cleaner->GetOutput()->GetNumberOfPoints();

		std::cout << "Original number of points " << Ntot << " after cleaning " << NCleaned << ". Removed " << Ntot-NCleaned << " points" << std::endl;

		// If there are not primitives associated with the points - cleaner eats them all
		if (NCleaned == 0 && Ntot > 0)
		{
			std::cerr << "Cleaner removed all points - we try with the originals instead" << std::endl;
		}

		vtkPolyDataPCANormals *PCANormals = vtkPolyDataPCANormals::New();
		if (NCleaned == 0)
			PCANormals->SetInputData(RawData);
		else
			PCANormals->SetInputData(cleaner->GetOutput());
		PCANormals->SetSearchRadius(m_Parms.NormalRadius);
		PCANormals->SetPointsPerNormals(m_Parms.NormalPointsPerNormal);
		PCANormals->SetSearchMode(m_Parms.NormalSearchMode);
		PCANormals->SetMax3EigenValue(m_Parms.Max3rdEigenval);
		PCANormals->SetMaxPlaneDistance(m_Parms.MaxPlaneDistance);
		//	PCANormals->SetScalarMode(VTK_SCALAR_MODE_PLANEDIST);
		PCANormals->Update();


		if (m_Parms.WriteLevel > 2)
		{
			vtkExtMisc::WritePDVTK(PCANormals->GetOutput(), outputnameNorms);
		}

		vtkPolyDataOrientNormalsByVoting *orientNormals = vtkPolyDataOrientNormalsByVoting::New();
		orientNormals->SetInputConnection(PCANormals->GetOutputPort());
		orientNormals->SetSearchRadius(m_Parms.CliqueNeighbourDistance);
		orientNormals->SetMinCliqueSize(m_Parms.MinCliqueSize);
		orientNormals->SetLargestCliqueOnly(m_Parms.LargestCliqueOnly);
		orientNormals->Update();

		// Simple check
		int NCheckPoints = orientNormals->GetOutput()->GetNumberOfPoints();

		if (NCheckPoints < 5)
		{
			std::cout << "Less than 5 points in oriented normals cloud" << std::endl;
			cleaner->Delete();
			PCANormals->Delete();
			orientNormals->Delete();
			return false;
		}
		NormalData = vtkPolyData::New();
		NormalData->DeepCopy(orientNormals->GetOutput());

		if (m_Parms.WriteLevel > 1)
		{
			vtkExtMisc::WritePDVTK(orientNormals->GetOutput(), outputnameONorms);
		}

		cleaner->Delete();
		PCANormals->Delete();
		orientNormals->Delete();

		CMultiTimer::MultiTimer().End("Estimating normals");
		return true;
	}
	else if (m_Parms.InputType == 1)
	{
		// treatments m1 and m3
		vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
		cleaner->SetInputData(RawData);
		cleaner->Update();

		int Ntot = RawData->GetNumberOfPoints();
		int NCleaned = cleaner->GetOutput()->GetNumberOfPoints();

		std::cout << "Original number of points " << Ntot << " after cleaning " << NCleaned << ". Removed " << Ntot-NCleaned << " points" << std::endl;
		// If there are not primitives associated with the points - cleaner eats them all
		if (NCleaned == 0 && Ntot > 0)
		{
			std::cerr << "Cleaner removed all points - we try with the originals instead" << std::endl;
		}

		vtkPolyDataPCANormals *PCANormals = vtkPolyDataPCANormals::New();
		if (NCleaned == 0)
			PCANormals->SetInputData(RawData);
		else
			PCANormals->SetInputConnection(cleaner->GetOutputPort());
		PCANormals->SetSearchRadius(m_Parms.NormalRadius);
		PCANormals->SetPointsPerNormals(m_Parms.NormalPointsPerNormal);
		PCANormals->SetSearchMode(m_Parms.NormalSearchMode);
		PCANormals->SetMax3EigenValue(m_Parms.Max3rdEigenval);
		PCANormals->SetMaxPlaneDistance(m_Parms.MaxPlaneDistance);
		//	PCANormals->SetScalarMode(VTK_SCALAR_MODE_PLANEDIST);
		PCANormals->Update();

		if (m_Parms.WriteLevel > 2)
		{
			vtkExtMisc::WritePDVTK(PCANormals->GetOutput(), outputnameNorms);
		}

		std::cout << "Using input data normals in the voting scheme" << std::endl;
		vtkPolyDataOrientNormalsByVoting *orientNormals = vtkPolyDataOrientNormalsByVoting::New();
		orientNormals->SetInputConnection(PCANormals->GetOutputPort());
		orientNormals->SetSearchRadius(m_Parms.CliqueNeighbourDistance);
		orientNormals->SetMinCliqueSize(m_Parms.MinCliqueSize);
		orientNormals->SetLargestCliqueOnly(m_Parms.LargestCliqueOnly);
		if (m_Parms.UseReferenceNormalsInVoting)
			orientNormals->SetReferencePolyData(RawData);
		orientNormals->Update();

		// Simple check
		int NCheckPoints = orientNormals->GetOutput()->GetNumberOfPoints();

		if (NCheckPoints < 5)
		{
			std::cout << "Less than 5 points in oriented normals cloud" << std::endl;
			cleaner->Delete();
			PCANormals->Delete();
			orientNormals->Delete();
			return false;
		}
		NormalData = vtkPolyData::New();
		NormalData->DeepCopy(orientNormals->GetOutput());

		if (m_Parms.WriteLevel > 1)
		{
			vtkExtMisc::WritePDVTK(orientNormals->GetOutput(), outputnameONorms);
		}

		cleaner->Delete();
		PCANormals->Delete();
		orientNormals->Delete();
		CMultiTimer::MultiTimer().End("Estimating normals");
		return true;
	}
	else if (m_Parms.InputType == 2 || m_Parms.InputType == 4 || m_Parms.InputType == 5)
	{
		// treatments (m1) and (m3)
		NormalData = vtkPolyData::New();

		if (m_Parms.EstimateNormals)
		{
			vtkPolyDataPCANormals *PCANormals = vtkPolyDataPCANormals::New();
			PCANormals->SetInputData(RawData);
			PCANormals->SetSearchRadius(m_Parms.NormalRadius);
			PCANormals->SetPointsPerNormals(m_Parms.NormalPointsPerNormal);
			PCANormals->SetSearchMode(m_Parms.NormalSearchMode);
			PCANormals->SetMax3EigenValue(m_Parms.Max3rdEigenval);
			PCANormals->SetMaxPlaneDistance(m_Parms.MaxPlaneDistance);
			//	PCANormals->SetScalarMode(VTK_SCALAR_MODE_PLANEDIST);
			PCANormals->Update();
			NormalData->DeepCopy(PCANormals->GetOutput());

			if (m_Parms.WriteLevel > 2)
			{
				vtkExtMisc::WritePDVTK(PCANormals->GetOutput(), outputnameNorms);
			}
			PCANormals->Delete();
		}
		else
		{
			NormalData->DeepCopy(RawData);
		}
		if (m_Parms.ConnectedComponentAnalysis)
		{
			std::cout << "Using input data normals in the voting scheme" << std::endl;
			vtkPolyDataOrientNormalsByVoting *orientNormals = vtkPolyDataOrientNormalsByVoting::New();
			orientNormals->SetInputData(NormalData);
			orientNormals->SetSearchRadius(m_Parms.CliqueNeighbourDistance);
			orientNormals->SetMinCliqueSize(m_Parms.MinCliqueSize);
			orientNormals->SetLargestCliqueOnly(m_Parms.LargestCliqueOnly);
			if (m_Parms.UseReferenceNormalsInVoting)
				orientNormals->SetReferencePolyData(RawData);
			orientNormals->Update();

			// Simple check
			int NCheckPoints = orientNormals->GetOutput()->GetNumberOfPoints();

			if (NCheckPoints < 5)
			{
				std::cout << "Less than 5 points in oriented normals cloud" << std::endl;
				orientNormals->Delete();
				NormalData->Delete();
				NormalData = NULL;
				return false;
			}
			NormalData->DeepCopy(orientNormals->GetOutput());
			orientNormals->Delete();
		}
		if (m_Parms.WriteLevel > 1)
		{
			vtkExtMisc::WritePDVTK(NormalData, outputnameONorms);
		}
		CMultiTimer::MultiTimer().End("Estimating normals");
		return true;
	}

	else if (m_Parms.InputType == 3)
	{
		NormalData = vtkPolyData::New();

		// treatments m4, (m1) and (m3)
		vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
		norms->SetInputData(RawData);
		norms->SplittingOff();
		norms->ConsistencyOn();
		norms->NonManifoldTraversalOff();
		norms->Update();

		if (m_Parms.EstimateNormals)
		{
			vtkPolyDataPCANormals *PCANormals = vtkPolyDataPCANormals::New();
			PCANormals->SetInputConnection(norms->GetOutputPort());
			PCANormals->SetSearchRadius(m_Parms.NormalRadius);
			PCANormals->SetPointsPerNormals(m_Parms.NormalPointsPerNormal);
			PCANormals->SetSearchMode(m_Parms.NormalSearchMode);
			PCANormals->SetMax3EigenValue(m_Parms.Max3rdEigenval);
			PCANormals->SetMaxPlaneDistance(m_Parms.MaxPlaneDistance);
			//	PCANormals->SetScalarMode(VTK_SCALAR_MODE_PLANEDIST);
			PCANormals->Update();
			NormalData->DeepCopy(PCANormals->GetOutput());

			if (m_Parms.WriteLevel > 2)
			{
				vtkExtMisc::WritePDVTK(PCANormals->GetOutput(), outputnameNorms);
			}
			PCANormals->Delete();
		}
		else
		{
			NormalData->DeepCopy(norms->GetOutput());
		}

		if (m_Parms.ConnectedComponentAnalysis)
		{
			std::cout << "Using input data normals in the voting scheme" << std::endl;
			vtkPolyDataOrientNormalsByVoting *orientNormals = vtkPolyDataOrientNormalsByVoting::New();
			orientNormals->SetInputData(NormalData);
			orientNormals->SetSearchRadius(m_Parms.CliqueNeighbourDistance);
			orientNormals->SetMinCliqueSize(m_Parms.MinCliqueSize);
			orientNormals->SetLargestCliqueOnly(m_Parms.LargestCliqueOnly);
			if (m_Parms.UseReferenceNormalsInVoting)
				orientNormals->SetReferencePolyData(RawData);
			orientNormals->Update();

			// Simple check
			int NCheckPoints = orientNormals->GetOutput()->GetNumberOfPoints();

			if (NCheckPoints < 5)
			{
				std::cout << "Less than 5 points in oriented normals cloud" << std::endl;
				orientNormals->Delete();
				NormalData->Delete();
				NormalData = NULL;
				norms->Delete();
				return false;
			}
			NormalData->DeepCopy(orientNormals->GetOutput());
			orientNormals->Delete();
		}
		if (m_Parms.WriteLevel > 1)
		{
			vtkExtMisc::WritePDVTK(NormalData, outputnameONorms);
		}

		norms->Delete();
		CMultiTimer::MultiTimer().End("Estimating normals");
		return true;
	}

	CMultiTimer::MultiTimer().End("Estimating normals");
	return true;
}

bool CMRFSurfaceReconstruction::ConvertRawScanData()
{
	std::string iname				= m_Parms.InputDir + m_Parms.InputName; 
	std::string outputname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_raw.vtk";

	if (RawData)
		return true;

	if (m_Parms.InputName != "" && file_exists(outputname.c_str()))
	{
		return true;
	}

	std::cout << "Converting raw scan" << std::endl;
	RawData = vtkExtMisc::MultiReadSurface(iname);
	if (!RawData)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return false;
	}

	// It is possible to rescale before processing
	if (m_Parms.ScaleInput)
	{
		std::cout << "Scaling input data" << std::endl;
		vtkPolyData *pdout  = vtkPolyData::New();
		double scale = 0;
		double CM[3];
		vtkExtMisc::NormaliseSize(RawData, pdout, 100, CM, scale);

		RawData->DeepCopy(pdout);
		pdout->Delete();
	}
	if (m_Parms.WriteLevel > 1)
	{
		vtkExtMisc::WritePDVTK(RawData, outputname);
	}

	return true;
}

void CMRFSurfaceReconstruction::TransformToOptimalBoundingBox()
{
	std::string BBname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_OrgCord.vtk";
	std::string iname1		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_OrgCord.vtk";
	std::string oname1		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string iname2		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_org_OrgCord.vtk";
	std::string oname2		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_org.vtk";
	std::string iname3		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_normals_OrgCord.vtk";
	std::string oname3		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_normals.vtk";
	std::string iname4		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_raw_OrgCord.vtk";
	std::string oname4		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_raw.vtk";

	if (m_Parms.InputName != "" && file_exists(oname1.c_str()))
	{
		return;
	}

	std::cout << "Transforming data to optimal bounding box" << std::endl;

	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(BBname.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << BBname << std::endl;
		return;
	}

	vtkPolyDataReader *readerI = vtkPolyDataReader::New();
	readerI->SetFileName(iname1.c_str());
	readerI->Update();
	if (!readerI)
	{
		std::cerr << "Could not read " << iname1 << std::endl;
		return;
	}

	double gridBounds[6];
	reader->GetOutput()->ComputeBounds();
	reader->GetOutput()->GetBounds(gridBounds);
	double Orgvolume = (gridBounds[1] - gridBounds[0]) * (gridBounds[3] - gridBounds[2]) * (gridBounds[5] - gridBounds[4]);

	std::cout << "Original BB volume " << Orgvolume << std::endl;

	vtkPrincipalAxisTransform *trans = vtkPrincipalAxisTransform::New();
	trans->SetSource(reader->GetOutput()->GetPoints());
	trans->SetDoScale(0);
	trans->Update();

	vtkTransform *tTrans = vtkTransform::New();
	tTrans->SetInput(trans);
	tTrans->Inverse();
	tTrans->Update();

	vtkTransformPolyDataFilter *filt = vtkTransformPolyDataFilter::New();
	filt->SetInputConnection(readerI->GetOutputPort());
	filt->SetTransform(tTrans);
	filt->Update();

	// Force it to point somewhat different
	vtkTransform *tTrans2 = vtkTransform::New();
	tTrans2->RotateY(180);
	tTrans2->RotateZ(90);
	tTrans2->Update();

	vtkTransformPolyDataFilter *filt2 = vtkTransformPolyDataFilter::New();
	filt2->SetInputConnection(filt->GetOutputPort());
	filt2->SetTransform(tTrans2);
	filt2->Update();

	filt2->GetOutput()->ComputeBounds();
	filt2->GetOutput()->GetBounds(gridBounds);
	double newvolume = (gridBounds[1] - gridBounds[0]) * (gridBounds[3] - gridBounds[2]) * (gridBounds[5] - gridBounds[4]);

	std::cout << "Transformed BB volume " << newvolume << std::endl;
	std::cout << "Volume reduced to " << newvolume / Orgvolume * 100 << " % of original volume" << std::endl;

	vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
	pWriter->SetInputConnection(filt2->GetOutputPort());
	pWriter->SetFileName(oname1.c_str());
	pWriter->Write();

	readerI->SetFileName(iname2.c_str());
	readerI->Update();
	filt2->Update();
	pWriter->SetFileName(oname2.c_str());
	pWriter->Update();
	pWriter->Write();


	readerI->SetFileName(iname3.c_str());
	pWriter->SetFileName(oname3.c_str());
	pWriter->Write();

	readerI->SetFileName(iname4.c_str());
	pWriter->SetFileName(oname4.c_str());
	pWriter->Write();

	trans->Delete();
	tTrans->Delete();
	filt->Delete();
	pWriter->Delete();
	reader->Delete();
	filt2->Delete();
	tTrans2->Delete();
}

void CMRFSurfaceReconstruction::CreateDensity()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string outputnameDens      = m_Parms.OutPutDir + m_Parms.InputNameShort + "_Density.vtk";

	if (m_Parms.UseLeavePatchOut)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}


	if (m_Parms.InputName != "" && file_exists(outputnameDens.c_str()))
	{
		return;
	}

	std::cout << "Creating density" << std::endl;


	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(iname.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	vtkPolyDataSamplingDensity *localDensity = vtkPolyDataSamplingDensity::New();
	localDensity->SetInputConnection(reader->GetOutputPort());
	localDensity->SetSearchRadius(m_Parms.NormalRadius);
	localDensity->Update();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(localDensity->GetOutputPort());
	writer->SetFileName(outputnameDens.c_str());
	writer->Write();

	writer->Delete();
	reader->Delete();

	localDensity->Delete();
}

void CMRFSurfaceReconstruction::ComputeDistanceMap()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string onameDistVol		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField.vti";
	std::string onameRefVol			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_ReferenceField.vti";
	std::string onameDistBounds		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.txt";
	std::string onameDistBoundsV	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.vtk";
	std::string onameDistDims		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Dimensions.txt";
	std::string onameDistCS			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_CellSize.txt";
	std::string onameVoxel			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_VoxelSize.vtk";
	std::string onameTimings		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Timings.txt";

	if (m_Parms.UseLeavePatchOut)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}

	if (m_Parms.InputName != "" && file_exists(onameDistVol.c_str()))
	{
		return;
	}

	vtkPolyDataReader *reader = vtkExtMisc::SafeReadPolyData(iname);
	if (!reader)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	vtkTimerLog *timer = vtkTimerLog::New();
	timer->StartTimer();

	std::cout << "Computing distance map" << std::endl;

	vtkOrientedPointSetDistanceFilter *signedDist = vtkOrientedPointSetDistanceFilter::New();
	signedDist->SetInputConnection(reader->GetOutputPort());
	signedDist->SetDistanceMode(m_Parms.DistanceMode);
	if (m_Parms.SampleSpace > 0)
		signedDist->SetSampleSpacing(m_Parms.SampleSpace);
	signedDist->SetSampleFactor(m_Parms.SampleFactor);
	signedDist->SetCreateReferenceVolume(m_Parms.UseLocalWeights);
	signedDist->SetNumberOfDistances(m_Parms.NumberOfDistances);
	signedDist->Update();

	timer->StopTimer();

	double RealSS = signedDist->GetSampleSpacing();

	double ElapsedTime = timer->GetElapsedTime();

	std::cout << "\nElapsedTime " << ElapsedTime << std::endl;

	std::cout << "Writing distance field" << std::endl;
	vtkXMLImageDataWriter* DSwriter = vtkXMLImageDataWriter::New();
	DSwriter->SetInputConnection(signedDist->GetOutputPort());
	DSwriter->SetFileName(onameDistVol.c_str());
	DSwriter->Write();
	DSwriter->Delete();

	if (m_Parms.UseLocalWeights)
	{
		std::cout << "Writing reference field" << std::endl;
		vtkXMLImageDataWriter* DSwriter2 = vtkXMLImageDataWriter::New();
		DSwriter2->SetInputData((vtkDataObject*)signedDist->GetReferenceVolume());
		DSwriter2->SetFileName(onameRefVol.c_str());
		DSwriter2->Write();
		DSwriter2->Delete();
	}

	double bounds[6];
	signedDist->GetOutput()->GetBounds(bounds);

	std::ofstream ost(onameDistBounds.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << onameDistBounds << std::endl;
	}
	else
	{
		ost << std::setprecision(10);
		for (int i = 0; i < 6; i++)
			ost << bounds[i] << " ";
	}

	std::ofstream ost2(onameDistCS.c_str());
	if (!ost2)
	{
		std::cerr << "Could not write to " << onameDistCS << std::endl;
	}
	else
	{
		ost2 << std::setprecision(10);
		ost2 << RealSS << std::endl;
	}

	int dims[3];
	signedDist->GetOutput()->GetDimensions(dims);

	std::ofstream ost3(onameDistDims.c_str());
	if (!ost3)
	{
		std::cerr << "Could not write to " << onameDistDims << std::endl;
	}
	else
	{
		ost3 << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
	}

	std::ofstream ost4(onameTimings.c_str());
	if (!ost4)
	{
		std::cerr << "Could not write to " << onameTimings << std::endl;
	}
	else
	{
		ost4 << ElapsedTime << std::endl;
	}

	vtkCubeSource *cube = vtkCubeSource::New();
	cube->SetBounds(bounds);

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(cube->GetOutputPort());
	writer->SetFileName(onameDistBoundsV.c_str());
	writer->Write();


	cube->SetCenter(0, 0, 0);
	cube->SetXLength(RealSS);
	cube->SetYLength(RealSS);
	cube->SetZLength(RealSS);
	cube->Update();

	writer->SetFileName(onameVoxel.c_str());
	writer->Write();


	writer->Delete();
	cube->Delete();


	signedDist->Delete();
	reader->Delete();
	timer->Delete();
}

void CMRFSurfaceReconstruction::MRFRegularisation()
{
	std::string inameDistVol		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField.vti";
	std::string inameRefVol			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_ReferenceField.vti";
	std::string inameDens			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Density.vtk";
	std::string onameWeigthVol		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_WeightField.vti";
	std::string onameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
	std::string onameTimings		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFRegularisation_Timings.txt";
	std::string onameValueTrack		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFRegularisation_ValueTrack.csv";

	if (m_Parms.InputName != "" && file_exists(onameMRFDist.c_str()))
	{
		return;
	}

	std::cout << "MRF based regularisation of the distance field" << std::endl;

	vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
	signedDist->SetFileName(inameDistVol.c_str());
	signedDist->Update();
	if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0)
	{
		std::cerr << "Could not read " << inameDistVol << std::endl;
		return;
	}

	vtkXMLImageDataReader *refVol = NULL;
	if (m_Parms.UseLocalWeights)
	{
		refVol = vtkXMLImageDataReader::New();
		refVol->SetFileName(inameRefVol.c_str());
		refVol->Update();
		if (!refVol)
		{
			std::cerr << "Could not read " << refVol << std::endl;
			return;
		}
	}

	vtkPolyDataReader *localDensity = vtkPolyDataReader::New();
	localDensity->SetFileName(inameDens.c_str());
	localDensity->Update();
	if (!localDensity)
	{
		std::cerr << "Could not read " << inameDens << std::endl;
		return;
	}

	vtkTimerLog *timer = vtkTimerLog::New();
	timer->StartTimer();

	//double space[3];
	//signedDist->GetOutput()->GetSpacing(space);
	//double RealSS = space[0];

	//double WeightValueHigh = std::max(m_Parms.LocalWeightMaxDist, RealSS*2);
	//std::cout << "Weight high value " << WeightValueHigh << std::endl;

	vtkDistanceFieldMRFRegularisation *filter = vtkDistanceFieldMRFRegularisation::New();
	filter->SetInputConnection(signedDist->GetOutputPort());
	if (m_Parms.UseLocalWeights)
		filter->SetReferenceVolume(refVol->GetOutput());
	else
		filter->SetReferenceVolume(NULL);
	filter->SetReferencePD(localDensity->GetOutput());
	filter->SetIterations(m_Parms.MaxIterations);
	filter->SetGlobalBeta(m_Parms.GlobalBeta);
	filter->SetWeightType(m_Parms.WeightMode);
	filter->SetPriorType(m_Parms.PriorType);
	filter->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
	filter->SetOptimisation(m_Parms.Optimisation);
	filter->Update();

	timer->StopTimer();

	double ElapsedTime = timer->GetElapsedTime();
	cout << "\nElapsedTime " << ElapsedTime << std::endl;

	std::cout << "Writing MRF distance field" << std::endl;
	vtkXMLImageDataWriter* DSwriter = vtkXMLImageDataWriter::New();
	DSwriter->SetInputData((vtkDataObject*)filter->GetOutput());
	DSwriter->SetFileName(onameMRFDist.c_str());
	DSwriter->Write();
	DSwriter->Delete();

	if (m_Parms.UseLocalWeights)
	{
		std::cout << "Writing weight field" << std::endl;
		vtkXMLImageDataWriter* DSwriter2 = vtkXMLImageDataWriter::New();
		DSwriter2->SetInputData((vtkDataObject*)filter->GetWeightVol());
		DSwriter2->SetFileName(onameWeigthVol.c_str());
		DSwriter2->Write();
		DSwriter2->Delete();
	}

	std::ofstream ost4(onameTimings.c_str());
	if (!ost4)
	{
		std::cerr << "Could not write to " << onameTimings << std::endl;
	}
	else
	{
		ost4 << ElapsedTime << std::endl;
	}

	std::vector<double> valtrack = filter->GetValueTrack();

	if (valtrack.size() > 0)
	{
		std::ofstream ost5(onameValueTrack.c_str());
		if (!ost5)
		{
			std::cerr << "Could not write to " << onameTimings << std::endl;
		}
		else
		{
			ost5 << std::setprecision(10);
			for (unsigned int i = 0; i < valtrack.size(); i++)
			{
				ost5 << valtrack[i] << std::endl;
			}
		}
	}


	filter->Delete();
	localDensity->Delete();
	signedDist->Delete();
	timer->Delete();
	if (refVol)
		refVol->Delete();
}

void CMRFSurfaceReconstruction::PolygoniseUpscaledDistanceMap()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string inameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
	std::string outputname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_Upscaled.vtk";

	if (m_Parms.InputName != "" && file_exists(outputname.c_str()))
	{
		return;
	}

	vtkPolyDataReader *reader =vtkExtMisc::SafeReadPolyData(iname);
	if (!reader)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	std::cout << "Recomputing distance map surface from distance map" << std::endl;

	vtkXMLImageDataReader *PreSD = vtkXMLImageDataReader::New();
	PreSD->SetFileName(inameMRFDist.c_str());
	PreSD->Update();
	if (!PreSD)
	{
		std::cerr << "Could not read " << inameMRFDist << std::endl;
		return;
	}

	vtkImageResample *resample = vtkImageResample::New();
	resample->SetInputConnection(PreSD->GetOutputPort());
	resample->SetAxisMagnificationFactor(0, 2);
	resample->SetAxisMagnificationFactor(1, 2);
	resample->SetAxisMagnificationFactor(2, 2);
	resample->Update();


	vtkOrientedPointSetDistanceFilter2 *signedDist = vtkOrientedPointSetDistanceFilter2::New();
	signedDist->SetInputConnection(reader->GetOutputPort());
	signedDist->SetInputDistances(resample->GetOutput());
	//	signedDist->SetDistanceMode(m_Parms.DistanceMode);
	//	if (m_Parms.SampleSpace > 0)
	//		signedDist->SetSampleSpacing(m_Parms.SampleSpace);
	//	signedDist->SetSampleFactor(m_Parms.SampleFactor);
	signedDist->SetCreateWeightVolume(true);
	signedDist->SetComputeMode(VTK_ORIENTEDDISTANCE_BAND);
	signedDist->SetNumberOfDistances(m_Parms.NumberOfDistances);
	signedDist->Update();

	vtkImplicitVolume *impvol = vtkImplicitVolume::New();
	impvol->SetVolume(signedDist->GetOutput());

	vtkPolyData *pd = vtkPolyData::New();
	std::string msg;
	double IsoValue = 0;
	if (PolygoniseAndRemesh(resample->GetOutput(), impvol, pd, m_Parms.Polygoniser, m_Parms.CellSizeFactor, m_Parms.Remesh, false, true, IsoValue, msg) != 0)
	{
		vtkExtMisc::WritePDVTK(pd, outputname);
	}
	std::cout << msg << std::endl;

	signedDist->Delete();
	pd->Delete();
	impvol->Delete();
	resample->Delete();
}


bool CMRFSurfaceReconstruction::PolygoniseDistanceMap()
{
	std::string iname               = m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string inameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
	std::string outputname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";

	if (MRFSurface)
		return true;

	if (m_Parms.InputName != "" && file_exists(outputname.c_str()))
	{
		return true;
	}

	std::cout << "Extracting isosurface from distance map" << std::endl;

	if (!MRFDistanceField)
	{
		vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
		signedDist->SetFileName(inameMRFDist.c_str());
		signedDist->Update();
		if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0)
		{
			std::cerr << "Could not read " << inameMRFDist << std::endl;
			return false;
		}
		MRFDistanceField = vtkImageData::New();
		MRFDistanceField->DeepCopy(signedDist->GetOutput());
		signedDist->Delete();
	}


	vtkImplicitVolume *impvol = vtkImplicitVolume::New();
	impvol->SetVolume(MRFDistanceField);

	MRFSurface = vtkPolyData::New();

	std::string msg;
	vtkPolyData *targetLenghtPd = NULL;
	if (m_Parms.RemeshToTargetEdgeLengths)
	{
		if (!NormalData)
		{
			NormalData = vtkExtMisc::MultiReadSurface(iname);
			if (!NormalData)
			{
				std::cerr << "Could not read " << iname << std::endl;
				return false;
			}
		}

		targetLenghtPd = NormalData;
	}
	double IsoValue = 0;
	if (m_Parms.SignedDistance == false && m_Parms.IsoValue == 0)
	{
		double spac[3];
		MRFDistanceField->GetSpacing(spac);
		double spacing = spac[0];
		IsoValue = spacing * 1.0;

		std::cout << "Estimated iso-value " << IsoValue << " from spacing " << spacing << std::endl;
	}
	else if (m_Parms.IsoValue != 0)
	{
		IsoValue = m_Parms.IsoValue;
		std::cout << "Iso-value " << IsoValue << " set from input parameter" << std::endl;
	}
	else
	{
		std::cout << "Iso-value " << IsoValue << " as default" << std::endl;
	}
	if (PolygoniseAndRemesh(MRFDistanceField, impvol, MRFSurface, m_Parms.Polygoniser, m_Parms.CellSizeFactor, m_Parms.Remesh, targetLenghtPd, true, IsoValue, msg) != 0)
	{
		if (m_Parms.WriteLevel > 1)
		{
			vtkExtMisc::WritePDVTK(MRFSurface, outputname);
		}
	}
	else
	{
		MRFSurface->Delete();
		MRFSurface = NULL;
		impvol->Delete();
		return false;
	}
	std::cout << msg << std::endl;

	impvol->Delete();
	return true;
}

void CMRFSurfaceReconstruction::PolygoniseRawDistanceMap()
{
	std::string inameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField.vti";
	std::string outputname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_RawSurface.vtk";

	if (m_Parms.InputName != "" && file_exists(outputname.c_str()))
	{
		return;
	}

	std::cout << "Extracting iso surface from raw distance map" << std::endl;

	vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
	signedDist->SetFileName(inameMRFDist.c_str());
	signedDist->Update();
	if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0)
	{
		std::cerr << "Could not read " << inameMRFDist << std::endl;
		return;
	}

	vtkImplicitVolume *impvol = vtkImplicitVolume::New();
	impvol->SetVolume(signedDist->GetOutput());

	vtkPolyData *pd = vtkPolyData::New();
	std::string msg;
	double IsoValue = 0;

	if (PolygoniseAndRemesh(signedDist->GetOutput(), impvol, pd, m_Parms.Polygoniser, m_Parms.CellSizeFactor,  m_Parms.Remesh, false, true, IsoValue, msg) != 0)
	{
		vtkExtMisc::WritePDVTK(pd, outputname);
	}
	std::cout << msg << std::endl;

	signedDist->Delete();
	pd->Delete();
	impvol->Delete();
}
void CMRFSurfaceReconstruction::CenterOfMass( vtkPoints* pd, double *CM )
{
	CM[0] = 0; CM[1] = 0; CM[2] = 0;

	const int N = pd->GetNumberOfPoints();

	if (!N) return;

	for (int i = 0; i < N; i++)
	{
		double pNi[3];
		pd->GetPoint(i, pNi);

		CM[0] += pNi[0];
		CM[1] += pNi[1];
		CM[2] += pNi[2];
	}
	CM[0] /= N; CM[1] /= N; CM[2] /= N;
}

bool CMRFSurfaceReconstruction::CropSurface()
{
	std::string pointname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string surfname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";
	std::string inameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
	//std::string inameDistBounds	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.txt";
	//std::string inameDistCS		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_CellSize.txt";
	std::string oname1			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_marked.vtk";
	std::string oname2			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
	std::string oname3			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_Aggresive_cropped.vtk";

	if (MRFSurfaceCropped)
		return true;

	if (m_Parms.UseLeavePatchOut)
	{
		pointname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		pointname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}

	if (m_Parms.InputName != "" && file_exists(oname1.c_str()))
	{
		return true;
	}
	std::cout << "Cropping surface" << std::endl;

	if (!NormalData)
	{
		NormalData = vtkExtMisc::MultiReadSurface(pointname);
		if (!NormalData)
		{
			std::cerr << "Could not read " << pointname << std::endl;
			return false;
		}
	}
	if (!MRFSurface)
	{
		MRFSurface = vtkExtMisc::MultiReadSurface(surfname);
		if (!MRFSurface)
		{
			std::cerr << "Could not read " << surfname << std::endl;
			return false;
		}
	}
	if (!MRFDistanceField)
	{
		vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
		signedDist->SetFileName(inameMRFDist.c_str());
		signedDist->Update();
		if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0)
		{
			std::cerr << "Could not read " << inameMRFDist << std::endl;
			return false;
		}
		MRFDistanceField = vtkImageData::New();
		MRFDistanceField->DeepCopy(signedDist->GetOutput());
		signedDist->Delete();
	}

	double bounds[6];
	MRFDistanceField->GetBounds(bounds);

	double spacing[3];
	MRFDistanceField->GetSpacing(spacing);
	double RealSS = spacing[0];

	CMultiTimer::MultiTimer().Start("Crop surface", 1);

	CMultiTimer::MultiTimer().Start("vtkPolyDataMarkCoveredCells", 1);
	vtkPolyDataMarkCoveredCells *marker = vtkPolyDataMarkCoveredCells::New();
	marker->SetInputData(MRFSurface);
	marker->SetSource(NormalData);
	marker->SetMarkRadius(m_Parms.NormalRadius);
	marker->Update();
	CMultiTimer::MultiTimer().End("vtkPolyDataMarkCoveredCells");

	CMultiTimer::MultiTimer().Start("vtkCropPolyData", 1);
	vtkCropPolyData *crop = vtkCropPolyData::New();
	crop->SetInputConnection(marker->GetOutputPort());
	crop->SetBoundary(bounds);
	crop->SetTolerance(RealSS);
	crop->Update();
	CMultiTimer::MultiTimer().End("vtkCropPolyData");

	if (m_Parms.WriteLevel > 2)
	{
		vtkExtMisc::WritePDVTK(crop->GetOutput(), oname1);
	}

	double threshold = 1.5;

	CMultiTimer::MultiTimer().Start("vtkClipPolyData", 1);
	vtkClipPolyData *clipper = vtkClipPolyData::New();
	clipper->SetInputConnection(crop->GetOutputPort());
	clipper->GenerateClipScalarsOff();
	clipper->GenerateClippedOutputOff();
	clipper->SetValue(threshold);
	clipper->InsideOutOn();
	clipper->Update();
	CMultiTimer::MultiTimer().End("vtkClipPolyData");

	if (m_Parms.WriteLevel > 0)
	{
		vtkExtMisc::WritePDVTK(clipper->GetOutput(), oname2);
	}

	MRFSurfaceCropped = vtkPolyData::New();
	MRFSurfaceCropped->DeepCopy(clipper->GetOutput());

	if (m_Parms.AggresiveCrop)
	{
		CMultiTimer::MultiTimer().Start("vtkClipPolyData", 1);
		// Keep only regions with value 1 (good support)
		double lowThreshold = 0.5;
		vtkClipPolyData *clipperLow = vtkClipPolyData::New();
		clipperLow->SetInputConnection(crop->GetOutputPort());
		clipperLow->GenerateClipScalarsOff();
		clipperLow->GenerateClippedOutputOff();
		clipperLow->SetValue(lowThreshold);
		clipperLow->InsideOutOff();
		clipperLow->Update();
		CMultiTimer::MultiTimer().End("vtkClipPolyData");

		CMultiTimer::MultiTimer().Start("vtkClipPolyData", 1);
		double highThreshold = 1.5;
		vtkClipPolyData *clipperHigh = vtkClipPolyData::New();
		clipperHigh->SetInputConnection(clipperLow->GetOutputPort());
		clipperHigh->GenerateClipScalarsOff();
		clipperHigh->GenerateClippedOutputOff();
		clipperHigh->SetValue(highThreshold);
		clipperHigh->InsideOutOn();
		clipperHigh->Update();
		CMultiTimer::MultiTimer().End("vtkClipPolyData");

		// Hack hack! 
		clipperHigh->GetOutput()->GetPointData()->SetScalars(NULL);

		if (m_Parms.WriteLevel > 0)
		{
			vtkExtMisc::WritePDVTK(clipperHigh->GetOutput(), oname3);
		}

		MRFSurfaceAggresivelyCropped = vtkPolyData::New();
		MRFSurfaceAggresivelyCropped->DeepCopy(clipperHigh->GetOutput());

		clipperLow->Delete();
		clipperHigh->Delete();
	}


	marker->Delete();
	crop->Delete();
	clipper->Delete();

	CMultiTimer::MultiTimer().End("Crop surface");

	return true;
}

vtkPolyData *CMRFSurfaceReconstruction::CropSurfaceExternal(vtkPolyData *RawData, vtkPolyData *Surface, const CParms& parms )
{
	std::cout << "Cropping surface" << std::endl;

	double bounds[6];
	Surface->GetBounds(bounds);

	vtkPolyDataMarkCoveredCells *marker = vtkPolyDataMarkCoveredCells::New();
	marker->SetInputData(Surface);
	marker->SetSource(RawData);
	marker->SetMarkRadius(parms.NormalRadius);
	marker->Update();

	vtkCropPolyData *crop = vtkCropPolyData::New();
	crop->SetInputConnection(marker->GetOutputPort());
	crop->SetBoundary(bounds);
	crop->SetTolerance(parms.MaxPlaneDistance * 1.5);
	crop->Update();

	vtkPolyData *pd = vtkPolyData::New();
	
	if (parms.AggresiveCrop)
	{
		// Keep only regions with value 1 (good support)
		double lowThreshold = 0.5;
		vtkClipPolyData *clipperLow = vtkClipPolyData::New();
		clipperLow->SetInputConnection(crop->GetOutputPort());
		clipperLow->GenerateClipScalarsOff();
		clipperLow->GenerateClippedOutputOff();
		clipperLow->SetValue(lowThreshold);
		clipperLow->InsideOutOff();
		clipperLow->Update();

		double highThreshold = 1.5;
		vtkClipPolyData *clipperHigh = vtkClipPolyData::New();
		clipperHigh->SetInputConnection(clipperLow->GetOutputPort());
		clipperHigh->GenerateClipScalarsOff();
		clipperHigh->GenerateClippedOutputOff();
		clipperHigh->SetValue(highThreshold);
		clipperHigh->InsideOutOn();
		clipperHigh->Update();

		// Hack hack! 
		clipperHigh->GetOutput()->GetPointData()->SetScalars(NULL);

		pd->DeepCopy(clipperHigh->GetOutput());
		clipperLow->Delete();
		clipperHigh->Delete();
	}
	else
	{
		double threshold = 1.5;
		vtkClipPolyData *clipper = vtkClipPolyData::New();
		clipper->SetInputConnection(crop->GetOutputPort());
		clipper->GenerateClipScalarsOff();
		clipper->GenerateClippedOutputOff();
		clipper->SetValue(threshold);
		clipper->InsideOutOn();
		clipper->Update();
		pd->DeepCopy(clipper->GetOutput());
		clipper->Delete();
	}

	marker->Delete();
	crop->Delete();

	pd->GetPointData()->SetScalars(NULL);
	return pd;
}


vtkPolyData * CMRFSurfaceReconstruction::DoApproximateColorSurface( vtkPolyData *RawData, vtkPolyData *Surface, const CParms& parms )
{
	vtkPolyDataDifference *diff = vtkPolyDataDifference::New();
	diff->SetInputData(Surface);
	diff->SetTargetData(RawData);
	diff->Update();

	vtkPolyData *pd = vtkPolyData::New();
	pd->DeepCopy(diff->GetOutput());

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetScalars());
	if (!scalars)
	{
		std::cerr << "No scalars!" << std::endl;
		pd->Delete();
		return NULL;
	}

	double mindist = parms.MaxPlaneDistance; // Default equals samplespacing
	double maxdist = parms.LocalWeightMaxDist;

	for (int i = 0; i < scalars->GetNumberOfTuples(); i++)
	{
		double val = scalars->GetValue(i);

		if (val < mindist)
			val = 0;
		if (val > maxdist)
			val = maxdist;
		val = val / maxdist;

		val = 1-val;

		scalars->SetValue(i, val);
	}

	scalars->Modified();
	pd->GetPointData()->SetScalars(scalars);

	diff->Delete();
	return pd;
}

void CMRFSurfaceReconstruction::ApproximateWeightColorSurface(int mode )
{
	std::string surfname		 = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";
	std::string iname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string oname			 = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_ApproximateWeights.vtk";

	// Do it on the cropped and remeshed one as well
	if (mode == 1)
	{
		surfname		 = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
		iname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
		oname			 = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_ApproximateWeights.vtk";
	}

	if (m_Parms.UseLeavePatchOut)
	{
		iname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		iname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}


	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	std::cout << "Approximate coloring surface based on weight" << std::endl;

	vtkPolyDataReader *RawData = vtkExtMisc::SafeReadPolyData(iname);
	if (!RawData)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	vtkPolyDataReader *Surface = vtkExtMisc::SafeReadPolyData(surfname);
	if (!Surface)
	{
		std::cerr << "Could not read " << surfname << std::endl;
		return;
	}

	vtkPolyData *pd = DoApproximateColorSurface(RawData->GetOutput(), Surface->GetOutput(), m_Parms);

	if (pd)
	{
		vtkExtMisc::WritePDVTK(pd, oname);
		pd->Delete();
	}

	Surface->Delete();
	RawData->Delete();
}

void CMRFSurfaceReconstruction::WeightColorSurface()
{
	std::string inameWeigthVol	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_WeightField.vti";
	std::string surfname		 = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";
	std::string oname			 = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_Weights.vtk";

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	std::cout << "Coloring surface based on weight" << std::endl;

	vtkXMLImageDataReader *weightVol = vtkXMLImageDataReader::New();
	weightVol->SetFileName(inameWeigthVol.c_str());
	weightVol->Update();
	if (!weightVol)
	{
		std::cerr << "Could not read " << inameWeigthVol << std::endl;
		return;
	}

	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(surfname.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << surfname << std::endl;
		return;
	}

	vtkPolyDataMarkUsingVolume *volmarker = vtkPolyDataMarkUsingVolume::New();
	volmarker->SetInputConnection(reader->GetOutputPort());
	volmarker->SetVolume(weightVol->GetOutput());
	volmarker->Update();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(volmarker->GetOutputPort());
	writer->SetFileName(oname.c_str());
	writer->Write();

	writer->Delete();
	volmarker->Delete();
	reader->Delete();
	weightVol->Delete();
}

bool CMRFSurfaceReconstruction::file_exists( char const* fn )
{
	struct stat fs;
	return stat(fn, &fs) == 0;
}

//void CMRFSurfaceReconstruction::OhtakeOptimiseMesh()
//{
//	std::string iname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
//	std::string inameMRFDist	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
//	std::string outputname		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_optimised.vtk";
//
//	if (file_exists(outputname.c_str()))
//	{
//		return;
//	}
//
//	std::cout << "Ohtake based optimising of mesh" << std::endl;
//
//	vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
//	signedDist->SetFileName(inameMRFDist.c_str());
//	signedDist->Update();
//	if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0))
//	{
//		std::cerr << "Could not read " << inameMRFDist << std::endl;
//		return;
//	}
//
//	vtkPolyDataReader *Surface = vtkPolyDataReader::New();
//	Surface->SetFileName(iname.c_str());
//	Surface->Update();
//	if (!Surface)
//	{
//		std::cerr << "Could not read " << iname << std::endl;
//		return;
//	}
//
//	vtkImplicitVolume *impvol = vtkImplicitVolume::New();
//	impvol->SetVolume(signedDist->GetOutputPort());
//
//	MeshData *mesh = new MeshData;
//	COhtakeUtils::PolyData2OhtakeMesh(Surface->GetOutputPort(), mesh);
//
//	Smoother *smoother = new Smoother;
//	smoother->setMeshData(mesh);
//
//	//	 Subdivider *subdivider = new Subdivider;
//	//	 subdivider->setMeshData(mesh);
//
//	//Simplifier *simplifier = new Simplifier;
//	//simplifier->setMeshData(mesh);
//
//
//	// 	 double subparm = 0.0005;
//	// 	 subdivider->subdivideCurvedFace(func, subparm);
//
//	// 	 mesh->computeFaceNormal();
//	// 	 mesh->computeNormal();
//
//	int iterations = 3;
//	double weight = 2;
//	double uniW = 0.001;
//	double ave = 0;
//
//	smoother->generateUniformCenter(impvol, iterations, weight);		
//	smoother->replaceOptimalPosition(500, uniW, ave);
//
//	mesh->computeFaceNormal();
//	mesh->computeNormal();
//
//	vtkPolyData *pd2 = vtkPolyData::New();
//	COhtakeUtils::OhtakeMesh2PolyData(mesh, pd2);
//
//	// Check gradients
//	//vtkPolyData *pd2 = vtkPolyData::New();
//	//vtkPoints *points = vtkPoints::New();
//	//vtkCellArray *lines = vtkCellArray::New();
//
//	//for (int i = 0; i < Surface->GetOutputPort()->GetNumberOfPoints(); i++)
//	//{
//	//	double p[3];
//	//	double n[3];
//	//	double pt[3];
//
//	//	Surface->GetOutputPort()->GetPoint(i, p);
//	//	
//	//	impvol->EvaluateGradient(p, n);
//
//	//	pt[0] = n[0] + p[0];
//	//	pt[1] = n[1] + p[1];
//	//	pt[2] = n[2] + p[2];
//
//	//	lines->InsertNextCell(2);
//	//	vtkIdType id = points->InsertNextPoint(p);
//	//	lines->InsertCellPoint(id);
//	//	id = points->InsertNextPoint(pt);
//	//	lines->InsertCellPoint(id);
//	//}
//
//	//pd2->SetPoints(points);
//	//points->Delete();
//	//pd2->SetLines(lines);
//	//lines->Delete();
//
//	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//	writer->SetInput(pd2);
//	writer->SetFileName(outputname.c_str());
//	writer->Write();
//
//	pd2->Delete();
//	impvol->Delete();
//	Surface->Delete();
//	signedDist->Delete();
//	//delete mesh;
//	//delete smoother;
//}

/*
%M
-0.758226351 -0.254590211 0.600230476 
-0.0188728304 0.928795505 0.370111774 
-0.651718203 0.269300452 -0.709042065 
%X 746.183145
%Y -253.248673
%Z 773.881994
*/

void CMRFSurfaceReconstruction::ReadCalibrationFiles()
{
	if (m_Parms.CalibrationDir == "")
		return;

	std::vector<std::string> CamNames;
	CamNames.push_back("1A");
	CamNames.push_back("1B");
	CamNames.push_back("1C");
	CamNames.push_back("2A");
	CamNames.push_back("2B");
	CamNames.push_back("2C");
	CamNames.push_back("3A");
	CamNames.push_back("3B");
	CamNames.push_back("3C");
	CamNames.push_back("4A");
	CamNames.push_back("4B");
	CamNames.push_back("4C");

	for (unsigned int i = 0; i < CamNames.size(); i++)
	{
		std::string name = CamNames[i];

		std::string iname = m_Parms.CalibrationDir + "calib_" + name + ".tka";
		std::string oname1b = m_Parms.OutPutDir + m_Parms.InputNameShort + "_Camera_" + name + ".vtk";

		C3dMDCamera cam;
		bool result = cam.Read(iname);
		if (result)
		{
			cam.ShowAsVTKFile(oname1b);
		}
	}
}

bool CMRFSurfaceReconstruction::AdaptiveParameters( vtkPolyData *pd, CParms &parms)
{
	const int NPoint = pd->GetNumberOfPoints();
	if (NPoint < 5)
	{
		std::cerr << "Less than 5 points" << std::endl;
		return false;
	}

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(pd);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	std::vector<double> Tdists;

	for (int i = 0; i < NPoint; i++)
	{
		double p[3];
		pd->GetPoint(i, p);

		vtkIdList *neighPts = vtkIdList::New();
		locator->FindClosestNPoints(2, p, neighPts);

		if (neighPts->GetNumberOfIds() == 2)
		{
			vtkIdType cid = neighPts->GetId(1);
			double cp[3];

			pd->GetPoint(cid, cp);


			double dist = sqrt(vtkMath::Distance2BetweenPoints(p, cp));
			if (dist > 0)
			{
				Tdists.push_back(dist);
			}
		}

		neighPts->Delete();
	}
	locator->Delete();

	int NPoints = pd->GetNumberOfPoints();
	double frac05 = 0;
	double frac95 = 0;
	double frac99 = 0;
	double median = 0;
	double mean = 0;
	double sdev = 0;
	double minScal = 0;
	double maxScal = 0;
	CGeneralUtils::MeanAndSdev(Tdists, mean, sdev);
	CGeneralUtils::MinMax(Tdists, minScal, maxScal);
	CGeneralUtils::Median(Tdists, 0.05, frac05);
	CGeneralUtils::Median(Tdists, 0.95, frac95);
	CGeneralUtils::Median(Tdists, 0.99, frac99);
	CGeneralUtils::Median(Tdists, 0.50, median);
	double RMS = CGeneralUtils::RMS(Tdists);

	std::cout << "Estimating parameters" << std::endl;
	parms.NormalRadius = parms.NormalRadiusFactor * mean;
	parms.MaxPlaneDistance = mean;
	parms.MinCliqueSize = NPoints * parms.MinCliquePercent;
	parms.CliqueNeighbourDistance = mean + 6.0 * sdev;
	parms.MinCellSize  = parms.SampleFactor * mean;
	parms.LocalWeightMaxDist = parms.kappa * mean;

	return true;
}


bool CMRFSurfaceReconstruction::InputDataStatistics()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_raw.vtk";
//	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_PointCloudStatistics.txt";
	std::string onameAdap			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_AdaptiveParameters.txt";

	if (!m_Parms.AdaptiveParameters && file_exists(oname.c_str()))
	{
		return true;
	}

	std::cout << "Computing point cloud statistics" << std::endl;

	if (!RawData)
	{
		RawData = vtkExtMisc::MultiReadSurface(iname);
		if (!RawData)
		{
			std::cerr << "Could not read " << iname << std::endl;
			return false;
		}
	}
	CMultiTimer::MultiTimer().Start("Estimating point cloud statistics 1", 1);
	std::cout << "Estimating point cloud statistics" << std::endl;

	const int NPoint = RawData->GetNumberOfPoints();
	std::vector<double> dists;

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(RawData);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	for (int i = 0; i < NPoint; i++)
	{
		double p[3];
		RawData->GetPoint(i, p);

		vtkIdList *neighPts = vtkIdList::New();
		locator->FindClosestNPoints(2, p, neighPts);

		if (neighPts->GetNumberOfIds() == 2)
		{
			vtkIdType cid = neighPts->GetId(1);
			double cp[3];

			RawData->GetPoint(cid, cp);

			double dist = sqrt(vtkMath::Distance2BetweenPoints(p, cp));
			if (dist > 0)
			{
				dists.push_back(dist);
			}
		}
		neighPts->Delete();
	}
	locator->Delete();
	CMultiTimer::MultiTimer().End("Estimating point cloud statistics 1");

	CMultiTimer::MultiTimer().Start("Estimating point cloud statistics 2", 1);
	int NPoints = RawData->GetNumberOfPoints();
	double frac05 = 0;
	double frac95 = 0;
	double frac99 = 0;
	double median = 0;
	double mean = 0;
	double sdev = 0;
	double minScal = 0;
	double maxScal = 0;
	CGeneralUtils::MeanAndSdev(dists, mean, sdev);
	CGeneralUtils::MinMax(dists, minScal, maxScal);
	CGeneralUtils::Median(dists, 0.05, frac05);
	CGeneralUtils::Median(dists, 0.95, frac95);
	CGeneralUtils::Median(dists, 0.99, frac99);
	CGeneralUtils::Median(dists, 0.50, median);
	double RMS = CGeneralUtils::RMS(dists);

	if (m_Parms.AdaptiveParameters)
	{
		std::cout << "Estimating adaptive parameters" << std::endl;
		m_Parms.NormalRadius = m_Parms.NormalRadiusFactor * mean;
		m_Parms.MaxPlaneDistance = mean;
		m_Parms.MinCliqueSize = NPoints * m_Parms.MinCliquePercent;
		m_Parms.CliqueNeighbourDistance = mean + 6.0 * sdev;
		m_Parms.MinCellSize  = m_Parms.SampleFactor * mean;
		m_Parms.LocalWeightMaxDist = m_Parms.kappa * mean;

		if (m_Parms.WriteLevel > 1)
		{
			std::ofstream ost(onameAdap.c_str());
			if (!ost)
			{
				std::cerr << "Could not write to " << onameAdap << std::endl;
			}
			else
			{
				ost << "NormalRadius "  << m_Parms.NormalRadius << std::endl;
				ost << "MaxPlaneDistance "  << m_Parms.MaxPlaneDistance << std::endl;
				ost << "MinCliqueSize "  << m_Parms.MinCliqueSize<< std::endl;
				ost << "CliqueNeighbourDistance "  << m_Parms.CliqueNeighbourDistance << std::endl;
				ost << "MinCellSize "  << m_Parms.MinCellSize << std::endl;
				ost << "LocalWeightMaxDist " << m_Parms.LocalWeightMaxDist << std::endl;
			}
		}
	}
	CMultiTimer::MultiTimer().End("Estimating point cloud statistics 2");

	if (m_Parms.WriteLevel > 1)
	{
		std::ofstream ost(oname.c_str());
		if (!ost)
		{
			std::cerr << "Could not write to " << oname << std::endl;
		}
		else
		{
			ost << "Point cloud neighbour point distances: " << std::endl;
			ost << " Min: "  << minScal << std::endl;
			ost << " Max: "  << maxScal << std::endl;
			ost << " Average: "  << mean << std::endl;
			ost << " Sdev: "  << sdev << std::endl;
			ost << " Median: "  << median << std::endl;
			ost << " RMS: "  << RMS << std::endl;
			ost << " 5% fractile: "  << frac05 << std::endl;
			ost << " 95% fractile: "  << frac95 << std::endl;
			ost << " 99% fractile: "  << frac99 << std::endl;
		}

	}

	return true;
}

void CMRFSurfaceReconstruction::ComputeAccuracy()
{
	std::string iname1				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string iname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_PointAccuracy.vtk";
	std::string oname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_PointAccuracy.txt";

	if (m_Parms.Remesh)
	{
		iname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
		oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_PointAccuracy.vtk";
		oname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_PointAccuracy.txt";
	}

	if (m_Parms.UseLeavePatchOut)
	{
		iname1		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}
	
	std::cout << "Computing accuracy" << std::endl;

	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(iname1.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << iname1 << std::endl;
		return;
	}

	vtkPolyDataReader *surf = vtkPolyDataReader::New();
	surf->SetFileName(iname2.c_str());
	surf->Update();
	if (!surf)
	{
		std::cerr << "Could not read " << iname2 << std::endl;
		return;
	}

	vtkPolyDataDifference *diff = vtkPolyDataDifference::New();
	diff->SetInputConnection(reader->GetOutputPort());
	diff->SetTargetData(surf->GetOutput());
//	diff->GetOutput();


	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(diff->GetOutputPort());
	writer->SetFileName(oname.c_str());
	writer->Write();


	std::ofstream fost(oname2.c_str());
	if (!fost)
	{
		std::cerr << "Could not write to " << oname2 << std::endl;
		return;
	}

	fost << vtkExtMisc::GetSurfaceValues(diff->GetOutput(),"\n"," ") << std::endl;

	writer->Delete();
	diff->Delete();
	surf->Delete();
	reader->Delete();
}

void CMRFSurfaceReconstruction::CreateCutSpheres()
{
	std::string iname 				= m_Parms.InputDir + m_Parms.InputNameShort + "_CutCoordinates.txt";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_CutSphere";

	std::ifstream fist(iname.c_str());
	if (!fist)
	{
		std::cerr << "Could not open " << iname << std::endl;
		return;
	}
	bool stop = false;
	int i = 1;
	do 
	{
		double p[3];
		double r = 0;
		
		fist >> p[0] >> p[1] >> p[2] >> r;
		if (fist.fail())
			stop = true;
		else
		{

			vtkSphereSource *sphere = vtkSphereSource::New();
			sphere->SetPhiResolution(50);
			sphere->SetThetaResolution(50);
			sphere->SetRadius(r);
			sphere->SetCenter(p);

			std::ostringstream ost;
			ost << oname << i << ".vtk";
//			std::cout << ost.str() << std::endl;

			vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
			writer->SetInputConnection(sphere->GetOutputPort());
			writer->SetFileName(ost.str().c_str());
			writer->Write();

			writer->Delete();
			sphere->Delete();
		}
		i++;
	} while (!stop);
}

void CMRFSurfaceReconstruction::CutWithSpheres()
{
	std::string sname 				= m_Parms.InputDir + m_Parms.InputNameShort + "_CutCoordinates.txt";
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_OrientedNormals_Cut.vtk";
	std::string oname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_OrientedNormals_Patch.vtk";

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	std::cout << "Cutting with spheres" << std::endl;

	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(iname.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	std::ifstream fist(sname.c_str());
	if (!fist)
	{
		std::cerr << "Could not open " << sname << std::endl;
		return;
	}

	vtkPolyData *pd = vtkPolyData::New();
	pd->DeepCopy(reader->GetOutput());

	vtkAppendPolyData *append = vtkAppendPolyData::New();

	bool stop = false;
	do 
	{
		double p[3];
		double r = 0;

		fist >> p[0] >> p[1] >> p[2] >> r;
		if (fist.fail())
			stop = true;
		else
		{
			vtkSphere *sphere = vtkSphere::New();
			sphere->SetRadius(r);
			sphere->SetCenter(p);
			sphere->Modified();

			vtkClipPolyData *clip = vtkClipPolyData::New();
			clip->SetInputData(pd);
			clip->SetClipFunction(sphere);
			clip->GenerateClippedOutputOn();
			clip->GenerateClipScalarsOff();
			clip->SetValue(0);
			clip->InsideOutOff();
			clip->Update();

			pd->DeepCopy(clip->GetOutput());

			// Arrghh! Why do append only works with temp pd2 ?
			vtkPolyData *pd2 = vtkPolyData::New();
			pd2->DeepCopy(clip->GetClippedOutput());
//			append->AddInput(clip->GetClippedOutput());
			append->AddInputData(pd2);
			append->Update();

			pd2->Delete();

			sphere->Delete();
			clip->Delete();
		}
	} while (!stop);

//	append->Modified();


	vtkRemoveUnusedPolyDataPoints *remover = vtkRemoveUnusedPolyDataPoints::New();
	remover->SetInputData(pd);
	remover->Update();

	vtkRemoveUnusedPolyDataPoints *remover2 = vtkRemoveUnusedPolyDataPoints::New();
	remover2->SetInputConnection(append->GetOutputPort());
//	remover2->SetInput(refFace);
	remover2->Update();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(remover->GetOutputPort());
	writer->SetFileName(oname.c_str());
	writer->Write();

	writer->SetInputConnection(remover2->GetOutputPort());
	writer->SetFileName(oname2.c_str());
	writer->Write();


	//	pd->Delete();
	remover->Delete();
	remover2->Delete();
	writer->Delete();
	append->Delete();
	pd->Delete();
//	refFace->Delete();
}	

//void CMRFSurfaceReconstruction::CutWithSpheres()
//{
//	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
//	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_OrientedNormals_Cut.vtk";
//	std::string oname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_OrientedNormals_Patch.vtk";
//
//	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
//	{
//		return;
//	}
//
//	std::cout << "Cutting with spheres" << std::endl;
//
//	vtkPolyDataReader *Surface = vtkPolyDataReader::New();
//	Surface->SetFileName(iname.c_str());
//	Surface->Update();
//	if (!Surface)
//	{
//		std::cerr << "Could not read " << iname << std::endl;
//		return;
//	}
//
//	double r1 = 25;
//	double r2 = 30;
//	double r3 = 30;
//	double c1[3];
//	c1[0] = 50; c1[1] = -25; c1[2] = 35; 
//	double c2[3];
//	c2[0] = -50; c2[1] = -90; c2[2] = 20; 
//	double c3[3];
//	c3[0] = -40; c3[1] = 90; c3[2] = 35; 
//
//	vtkSphere *sphere1 = vtkSphere::New();
//	sphere1->SetRadius(r1);
//	sphere1->SetCenter(c1);
//	sphere1->Modified();
//
//	vtkSphere *sphere2 = vtkSphere::New();
//	sphere2->SetRadius(r2);
//	sphere2->SetCenter(c2);
//	sphere2->Modified();
//
//	//vtkClipPolyData *clipperLow = vtkClipPolyData::New();
//	//clipperLow->SetInput(ptransSource->GetOutputPort());
//	//clipperLow->SetClipFunction(cplane);
//	//clipperLow->GenerateClipScalarsOff();
//	//clipperLow->GenerateClippedOutputOff();
//	//clipperLow->SetValue(0);
//	//clipperLow->InsideOutOff();
//	//clipperLow->Update();
//
//	vtkClipPolyData *clip = vtkClipPolyData::New();
//	clip->SetInput(Surface->GetOutputPort());
//	clip->SetClipFunction(sphere1);
//	clip->GenerateClippedOutputOn();
//	clip->GenerateClipScalarsOff();
//	clip->SetValue(0);
//	clip->InsideOutOff();
//	clip->Update();
//
//// 	vtkPolyData *pd = vtkPolyData::New();
//// 	pd->DeepCopy(clip->GetOutputPort());
//
//	//vtkRemoveUnusedPolyDataPoints *remover = vtkRemoveUnusedPolyDataPoints::New();
//	//remover->SetInput(clip->GetOutputPort());
//	//remover->Update();
//
//	vtkAppendPolyData *append = vtkAppendPolyData::New();
//	
//	//remover->SetInput(clip->GetClippedOutput());
//	//remover->Update();
//
//	append->AddInput(clip->GetClippedOutput());
//
//	vtkClipPolyData *clip2 = vtkClipPolyData::New();
//	clip2->SetInput(clip->GetOutputPort());
//	clip2->SetClipFunction(sphere2);
//	clip2->GenerateClippedOutputOn();
//	clip2->GenerateClipScalarsOff();
//	clip2->SetValue(0);
//	clip2->InsideOutOff();
//	clip2->Update();
//
//	append->AddInput(clip2->GetClippedOutput());
//
//	vtkRemoveUnusedPolyDataPoints *remover = vtkRemoveUnusedPolyDataPoints::New();
//	remover->SetInput(clip2->GetOutputPort());
//	remover->Update();
//
//	vtkRemoveUnusedPolyDataPoints *remover2 = vtkRemoveUnusedPolyDataPoints::New();
//	remover2->SetInput(append->GetOutputPort());
//	remover2->Update();
//
//
//	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//	//	writer->SetInput(clip->GetOutputPort());
//	//	writer->SetInput(pd);
//	writer->SetInput(remover->GetOutputPort());
//	writer->SetFileName(oname.c_str());
//	writer->Write();
//
//	writer->SetInput(remover2->GetOutputPort());
//	writer->SetFileName(oname2.c_str());
//	writer->Write();
//
//
//
////	pd->Delete();
//	remover->Delete();
//	remover2->Delete();
//	writer->Delete();
//	clip->Delete();
//	clip2->Delete();
//	sphere1->Delete();
//	sphere2->Delete();
//	append->Delete();
//}	

void CMRFSurfaceReconstruction::ComputePatchAccuracy()
{
	std::string iname1				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_patch.vtk";
	std::string iname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_Patch_PointAccuracy.vtk";
	std::string oname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_Patch_PointAccuracy.txt";

	if (m_Parms.Remesh)
	{
		iname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
		oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_Patch_PointAccuracy.vtk";
		oname2				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_Patch_PointAccuracy.txt";
	}

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	std::cout << "Computing patch accuracy" << std::endl;

	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(iname1.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << iname1 << std::endl;
		return;
	}

	vtkPolyDataReader *surf = vtkPolyDataReader::New();
	surf->SetFileName(iname2.c_str());
	surf->Update();
	if (!surf)
	{
		std::cerr << "Could not read " << iname2 << std::endl;
		return;
	}

	vtkPolyDataDifference *diff = vtkPolyDataDifference::New();
	diff->SetInputConnection(reader->GetOutputPort());
	diff->SetTargetData(surf->GetOutput());
//	diff->GetOutput();


	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(diff->GetOutputPort());
	writer->SetFileName(oname.c_str());
	writer->Write();


	std::ofstream fost(oname2.c_str());
	if (!fost)
	{
		std::cerr << "Could not write to " << oname2 << std::endl;
		return;
	}

	fost << vtkExtMisc::GetSurfaceValues(diff->GetOutput(),"\n"," ") << std::endl;

	writer->Delete();
	diff->Delete();
	surf->Delete();
	reader->Delete();
}

void CMRFSurfaceReconstruction::AddNoiseToCloud()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_OrientedNormals_Noise.vtk";

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	std::cout << "Adding noise" << std::endl;

	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(iname.c_str());
	reader->Update();
	if (!reader)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}


	vtkPolyData *pd = vtkPolyData::New();
	pd->DeepCopy(reader->GetOutput());


	// Type of noise
	// 0 outlier noise, 1 Gaussian on all points1
	if (m_Parms.NoiseNature == 1)
	{
		std::cout << "NoiseType Gaussian" << std::endl;
		double noiseSdev = 1.0;
		for (int i = 0; i < pd->GetNumberOfPoints(); i++)
		{
			double p[3];
			pd->GetPoint(i, p);

			p[0] = p[0] + CGeneralUtils::GaussRandomNumber() * noiseSdev;
			p[1] = p[1] + CGeneralUtils::GaussRandomNumber() * noiseSdev;
			p[2] = p[2] + CGeneralUtils::GaussRandomNumber() * noiseSdev;

			pd->GetPoints()->SetPoint(i, p);
		}
	}
	else
	{
		std::cout << "NoiseType Outlier" << std::endl;
		// probability of an outlier
		double OutlierProbability = 5.0;
		double distMax = 3;

		int setnum = 0;

		for (int i = 0; i < pd->GetNumberOfPoints(); i++)
		{
			double num = vtkMath::Random(0, 100);
			if (num < OutlierProbability)
			{
				double p[3];
				pd->GetPoint(i, p);

				p[0] = p[0] + vtkMath::Random(-distMax, distMax);
				p[1] = p[1] + vtkMath::Random(-distMax, distMax);
				p[2] = p[2] + vtkMath::Random(-distMax, distMax);

				pd->GetPoints()->SetPoint(i, p);
				setnum++;
			}
		}
		std::cout << "Moved " << setnum << " out of " << pd->GetNumberOfPoints() << " points" << std::endl;
	}

	pd->Modified();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputData(pd);
	writer->SetFileName(oname.c_str());
	writer->Write();

	writer->Delete();
	reader->Delete();
}


void CMRFSurfaceReconstruction::SparseCholeskyMRFRegularisation()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField_Initial.vti";

	std::string onameDistBounds		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Initial_Bounds.txt";
	std::string onameDistBoundsV	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Initial_Bounds.vtk";
	std::string onameDistDims		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Initial_Dimensions.txt";
	std::string onameDistCS			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Initial_CellSize.txt";
	std::string onameVoxel			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Initial_VoxelSize.vtk";

	if (m_Parms.UseLeavePatchOut)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	vtkPolyDataReader *surface = vtkPolyDataReader::New();
	surface->SetFileName(iname.c_str());
	surface->Update();
	if (!surface)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	std::cout << "Computing initial distance map" << std::endl;

	int minSideVoxels = 50;
	int MaxNumberOfVoxels = 200000;

	vtkOrientedPointSetDistanceFilter2 *SDF = vtkOrientedPointSetDistanceFilter2::New();
	SDF->SetInputConnection(surface->GetOutputPort());
	SDF->SetDistanceMode(m_Parms.DistanceMode);
	SDF->SetCreateWeightVolume(true);
	SDF->SetMinSideVoxels(minSideVoxels);
	SDF->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
	SDF->SetNumberOfDistances(m_Parms.NumberOfDistances);
	SDF->SetPadVoxels(m_Parms.PadVoxels);
	SDF->SetComputeMode(VTK_ORIENTEDDISTANCE_COMBINED_FULLBAND);
	SDF->Update();

	vtkImageData *WeigthVolume = SDF->GetWeightVolume();

	std::cout << "MRF Regularisation" << std::endl;

	//double space[3];
	//SDF->GetOutput()->GetSpacing(space);
	//double RealSS = space[0];

	//double WeightValueHigh = std::max(m_Parms.LocalWeightMaxDist, RealSS*2);
	//std::cout << "Weight high value " << WeightValueHigh << std::endl;


	vtkDistanceFieldMRFRegularisation *filter = vtkDistanceFieldMRFRegularisation::New();
	filter->SetInputConnection(SDF->GetOutputPort());
	filter->SetWeightVol(WeigthVolume);
	//	filter->SetIterations(m_Parms.MaxIterations);
	filter->SetIterations(1000);
	filter->SetGlobalBeta(m_Parms.GlobalBeta);
	filter->SetWeightType(m_Parms.WeightMode);
	filter->SetPriorType(m_Parms.PriorType);
	filter->SetOptimisation(m_Parms.Optimisation);
	filter->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
	filter->Update();

	vtkDistanceFieldMRFRegularisation *resultMRF = filter;

	std::cout << "Writing MRF distance field" << std::endl;
	vtkXMLImageDataWriter* DSwriter = vtkXMLImageDataWriter::New();
	DSwriter->SetInputData((vtkDataObject*)resultMRF->GetOutput());
	DSwriter->SetFileName(oname.c_str());
	DSwriter->Write();

	double bounds[6];
	resultMRF->GetOutput()->GetBounds(bounds);

	double spacing[3];
	resultMRF->GetOutput()->GetSpacing(spacing);
	double RealSS = spacing[0];

	std::ofstream ost(onameDistBounds.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << onameDistBounds << std::endl;
	}
	else
	{
		ost << std::setprecision(10);
		for (int i = 0; i < 6; i++)
			ost << bounds[i] << " ";
	}

	std::ofstream ost2(onameDistCS.c_str());
	if (!ost2)
	{
		std::cerr << "Could not write to " << onameDistCS << std::endl;
	}
	else
	{
		ost2 << std::setprecision(10);
		ost2 << RealSS << std::endl;
	}

	int dims[3];
	resultMRF->GetOutput()->GetDimensions(dims);

	std::ofstream ost3(onameDistDims.c_str());
	if (!ost3)
	{
		std::cerr << "Could not write to " << onameDistDims << std::endl;
	}
	else
	{
		ost3 << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
	}

	vtkCubeSource *cube = vtkCubeSource::New();
	cube->SetBounds(bounds);

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(cube->GetOutputPort());
	writer->SetFileName(onameDistBoundsV.c_str());
	writer->Write();

	cube->SetCenter(0, 0, 0);
	cube->SetXLength(RealSS);
	cube->SetYLength(RealSS);
	cube->SetZLength(RealSS);
	cube->Update();

	writer->SetFileName(onameVoxel.c_str());
	writer->Write();

	writer->Delete();
	cube->Delete();

	SDF->Delete();
	DSwriter->Delete();

	filter->Delete();
	surface->Delete();
}

void CMRFSurfaceReconstruction::MultiLevelMRFRefinement()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string inameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField_Initial.vti";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";

	std::string onameDistBounds		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.txt";
	std::string onameDistBoundsV	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.vtk";
	std::string onameDistDims		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Dimensions.txt";
	std::string onameDistCS			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_CellSize.txt";
	std::string onameVoxel			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_VoxelSize.vtk";

	if (m_Parms.UseLeavePatchOut)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}

	vtkPolyDataReader *surface = vtkPolyDataReader::New();
	surface->SetFileName(iname.c_str());
	surface->Update();
	if (!surface)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return;
	}

	vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
	signedDist->SetFileName(inameMRFDist.c_str());
	signedDist->Update();
	if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0)
	{
		std::cerr << "Could not read " << inameMRFDist << std::endl;
		return;
	}

	std::cout << "Multi-level refinement of distance map" << std::endl;


	int maxLevels = 20;

	std::vector<vtkImageResample*> Resamplers(maxLevels, NULL);
	std::vector<vtkOrientedPointSetDistanceFilter2*> SD(maxLevels, NULL);
	std::vector<vtkDistanceFieldMRFRegularisation*> filters(maxLevels, NULL);

	vtkImageData* curVolume = signedDist->GetOutput();

	int iterations = 50;
	int MaxVolumeSize = m_Parms.MaxVolumeSize;
	int id = 0;
	bool stop = false;
	do 
	{
		std::cout << "*****---- Iteration " << id+1 << " ----*****" << std::endl;
		DWORDLONG freeDoubles = CMemoryUtils::GetAvailableMemoryForDouble();
		std::cout << "System claims space for " << freeDoubles << " doubles" << std::endl;
		std::cout << "Upsampling field" << std::endl;

		Resamplers[id] = vtkImageResample::New();

		Resamplers[id]->SetInputData(curVolume);
		Resamplers[id]->SetAxisMagnificationFactor(0, 2);
		Resamplers[id]->SetAxisMagnificationFactor(1, 2);
		Resamplers[id]->SetAxisMagnificationFactor(2, 2);
		Resamplers[id]->Update();

		int resDim[3];
		Resamplers[id]->GetOutput()->GetDimensions(resDim);
		int resSize = resDim[0] * resDim[1] * resDim[2]; 
		if (resSize > MaxVolumeSize)
		{
			std::cout << "Size limit reached with (" << resDim[0] << ", " << resDim[1] << ", " << resDim[2] << ") " << resSize << std::endl;
			stop = true;
		}
		else
		{
			std::cout << "Size of current volume (" << resDim[0] << ", " << resDim[1] << ", " << resDim[2] << ") " << resSize << std::endl;
			std::cout << "Calculating banded distance" << std::endl;

			SD[id] = vtkOrientedPointSetDistanceFilter2::New();
			SD[id]->SetInputConnection(surface->GetOutputPort());
			SD[id]->SetInputDistances(Resamplers[id]->GetOutput());
			SD[id]->SetDistanceMode(m_Parms.DistanceMode);
			SD[id]->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
			SD[id]->SetCreateWeightVolume(true);
			SD[id]->SetComputeMode(VTK_ORIENTEDDISTANCE_BAND);
			SD[id]->SetNumberOfDistances(m_Parms.NumberOfDistances);
			SD[id]->SetPadVoxels(m_Parms.PadVoxels);
			SD[id]->Update();

			vtkImageData *WeigthVolume = SD[id]->GetWeightVolume();

			std::cout << "Regularisation" << std::endl;
			filters[id] = vtkDistanceFieldMRFRegularisation::New();
			filters[id]->SetInputConnection(SD[id]->GetOutputPort());
			filters[id]->SetWeightVol(WeigthVolume);
			filters[id]->SetGlobalBeta(m_Parms.GlobalBeta);
			filters[id]->SetWeightType(m_Parms.WeightMode);
			filters[id]->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
			filters[id]->SetPriorType(m_Parms.PriorType);
			filters[id]->SetOptimisation(0);
			//	filter2->SetIterations(m_Parms.MaxIterations);
			filters[id]->SetIterations(iterations);
			filters[id]->SetMinRMS(0.01);
			//filters[id]->SetCGTolerance(1);
			filters[id]->Update();

			curVolume = filters[id]->GetOutput();
		}

		iterations /= 4;

		id++;
		if (id >= maxLevels)
			stop = true;

	} while (!stop);

	std::cout << "Writing MRF distance field" << std::endl;
	vtkXMLImageDataWriter* DSwriter = vtkXMLImageDataWriter::New();
	DSwriter->SetInputData((vtkDataObject*)curVolume);
	DSwriter->SetFileName(oname.c_str());
	DSwriter->Write();
	DSwriter->Delete();

	double bounds[6];
	curVolume->GetBounds(bounds);

	double spacing[3];
	curVolume->GetSpacing(spacing);
	double RealSS = spacing[0];

	std::ofstream ost(onameDistBounds.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << onameDistBounds << std::endl;
	}
	else
	{
		ost << std::setprecision(10);
		for (int i = 0; i < 6; i++)
			ost << bounds[i] << " ";
	}

	std::ofstream ost2(onameDistCS.c_str());
	if (!ost2)
	{
		std::cerr << "Could not write to " << onameDistCS << std::endl;
	}
	else
	{
		ost2 << std::setprecision(10);
		ost2 << RealSS << std::endl;
	}

	int dims[3];
	curVolume->GetDimensions(dims);

	std::ofstream ost3(onameDistDims.c_str());
	if (!ost3)
	{
		std::cerr << "Could not write to " << onameDistDims << std::endl;
	}
	else
	{
		ost3 << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
	}

	vtkCubeSource *cube = vtkCubeSource::New();
	cube->SetBounds(bounds);

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(cube->GetOutputPort());
	writer->SetFileName(onameDistBoundsV.c_str());
	writer->Write();


	cube->SetCenter(0, 0, 0);
	cube->SetXLength(RealSS);
	cube->SetYLength(RealSS);
	cube->SetZLength(RealSS);
	cube->Update();

	writer->SetFileName(onameVoxel.c_str());
	writer->Write();

	writer->Delete();
	cube->Delete();
	surface->Delete();
	signedDist->Delete();

	for (int i = 0; i < maxLevels; i++)
	{
		if (Resamplers[i])
			Resamplers[i]->Delete();
		if (SD[i])
			SD[i]->Delete();
		if (filters[i])
			filters[i]->Delete();
	}
}


void CMRFSurfaceReconstruction::WriteDiagnosisFile( double SplineEnergy, int iteration )
{
	std::ostringstream name;
	name << m_Parms.OutPutDir << m_Parms.InputNameShort << "_Diagnostics_it"<< iteration << ".txt";

	std::ofstream fost(name.str().c_str());

	fost << "Final SplineEnergy: " << SplineEnergy << std::endl;
}

void CMRFSurfaceReconstruction::WriteBoundsFile( vtkImageData *dat, int iteration )
{
	std::ostringstream name;
	name << m_Parms.OutPutDir << m_Parms.InputNameShort << "_DistanceField_Bounds_it"<< iteration << ".vtk";

	double bounds[6];
	dat->GetBounds(bounds);

	vtkCubeSource *cube = vtkCubeSource::New();
	cube->SetBounds(bounds);

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(cube->GetOutputPort());
	writer->SetFileName(name.str().c_str());
	writer->Write();

	cube->Delete();
	writer->Delete();

	int dims[3];
	dat->GetDimensions(dims);
	name.str("");
	name << m_Parms.OutPutDir << m_Parms.InputNameShort << "_DistanceField_Dimensions_it"<< iteration << ".txt";

	std::ofstream ost(name.str().c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << name.str() << std::endl;
	}
	else
	{
		ost << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
	}

	name.str("");
	name << m_Parms.OutPutDir << m_Parms.InputNameShort << "_DistanceField_Bounds_it"<< iteration << ".txt";

	std::ofstream ost2(name.str().c_str());
	if (!ost2)
	{
		std::cerr << "Could not write to " << name.str() << std::endl;
	}
	else
	{
		ost2 << std::setprecision(10);
		for (int i = 0; i < 6; i++)
			ost2 << bounds[i] << " ";
	}

	double spacing[3];
	dat->GetSpacing(spacing);
	double RealSS = spacing[0];

	name.str("");
	name << m_Parms.OutPutDir << m_Parms.InputNameShort << "_DistanceField_Spacing_it"<< iteration << ".txt";

	std::ofstream ost3(name.str().c_str());
	if (!ost3)
	{
		std::cerr << "Could not write to " << name.str() << std::endl;
	}
	else
	{
		ost3 << std::setprecision(10);
		ost3 << RealSS << std::endl;
	}


}

void CMRFSurfaceReconstruction::WriteVoxelSizeFile( vtkImageData *dat, int iteration )
{
	std::ostringstream name;
	name << m_Parms.OutPutDir << m_Parms.InputNameShort << "_DistanceField_VoxelSize_it"<< iteration << ".vtk";

	double spacing[3];
	dat->GetSpacing(spacing);
	double RealSS = spacing[0];

	vtkCubeSource *cube = vtkCubeSource::New();
	cube->SetCenter(0, 0, 0);
	cube->SetXLength(RealSS);
	cube->SetYLength(RealSS);
	cube->SetZLength(RealSS);
	cube->Update();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputConnection(cube->GetOutputPort());
	writer->SetFileName(name.str().c_str());
	writer->Write();

	cube->Delete();
	writer->Delete();
}

int CMRFSurfaceReconstruction::EstimateNewVolumeDimensionsWithoutShrink( vtkImageData *vol, vtkPolyData *points, int *extent, double &sampleSpacing, double *origin, int maxSize, double MinCellSize )
{
	// Get distance between voxels in org volume (cellsize)
	double spacing[3];
	vol->GetSpacing(spacing);
	double RealSS = spacing[0];

	//  half the samplespacing
	double ss = RealSS / 2;

	// Get dimensions of org volume
	double VolBounds[6];
	vol->GetBounds(VolBounds);

	double Xs = VolBounds[1] - VolBounds[0];
	double Ys = VolBounds[3] - VolBounds[2];
	double Zs = VolBounds[5] - VolBounds[4];

	int sx = (int)(Xs/ss);
	int sy = (int)(Ys/ss);
	int sz = (int)(Zs/ss);

	int size = sx * sy * sz;

	sampleSpacing = ss;
	double Xmin = VolBounds[0];
	double Ymin = VolBounds[2];
	double Zmin = VolBounds[4];

	origin[0] = Xmin;
	origin[1] = Ymin;
	origin[2] = Zmin;

	extent[0] = 0;
	extent[1] = sx;
	extent[2] = 0;
	extent[3] = sy;
	extent[4] = 0;
	extent[5] = sz;

	return size;
}

void CMRFSurfaceReconstruction::CheckVolumeIfIsoSurfacePresentAtEdges(vtkImageData *vol, vtkPolyData *points, int *extent, double &sampleSpacing, double *origin, int maxSize, double MinCellSize)
{
	// First try slice in top 
	// We compute real coordinate of ymax (suggested new volume size)
	double y = extent[3] * sampleSpacing + origin[1];

	double vbounds[6];
	vol->GetBounds(vbounds);
	double spacing[3];
	vol->GetSpacing(spacing);

	int dims[3] = { 0,0,0 };
	vol->GetDimensions(dims);

	// Now compute slice id
	int yslice = (y - vbounds[2]) / spacing[0];

	if (yslice < 0 || yslice > dims[1])
	{
		std::cerr << "Something wrong with y coord in CheckVolumeIfIsoSurfacePresentAtEdges" << std::endl;
		return;
	}

	vtkDoubleArray *values =
		vtkDoubleArray::SafeDownCast(vol->GetPointData()->GetScalars());

	bool negativeVal = false;
	bool positiveVal = false;
	int yoff = yslice * dims[0];
	for (int x = 0; x < dims[0]; x++)
	{
		for (int z = 0; z < dims[2]; z++)
		{
			int offset = z * dims[0] * dims[1] + yoff + x;
			double val = values->GetValue(offset);
			if (val > 0)
				positiveVal = true;
			if (val < 0)
				negativeVal = true;
			if (positiveVal && negativeVal)
				break;
		}
		if (positiveVal && negativeVal)
			break;
	}

	if (positiveVal && negativeVal)
		std::cout << "iso-surface found at y max" << std::endl;
	else
		std::cout << "NO iso-surface found at y max" << std::endl;
}

bool CMRFSurfaceReconstruction::CheckYSliceIsoSurface(vtkImageData *vol, int y)
{
	int dims[3] = { 0,0,0 };
	vol->GetDimensions(dims);

	vtkDoubleArray *values =
		vtkDoubleArray::SafeDownCast(vol->GetPointData()->GetScalars());

	bool negativeVal = false;
	bool positiveVal = false;
	int yoff = y * dims[0];
	for (int x = 0; x < dims[0]; x++)
	{
		for (int z = 0; z < dims[2]; z++)
		{
			int offset = z * dims[0] * dims[1] + yoff + x;
			double val = values->GetValue(offset);
			if (val > 0)
				positiveVal = true;
			if (val < 0)
				negativeVal = true;
			if (positiveVal && negativeVal)
				break;
		}
		if (positiveVal && negativeVal)
			break;
	}

	bool IsoSurf = positiveVal && negativeVal;

	return IsoSurf;
}

bool CMRFSurfaceReconstruction::CheckXSliceIsoSurface(vtkImageData *vol, int x)
{
	int dims[3] = { 0,0,0 };
	vol->GetDimensions(dims);

	vtkDoubleArray *values =
		vtkDoubleArray::SafeDownCast(vol->GetPointData()->GetScalars());

	bool negativeVal = false;
	bool positiveVal = false;
	for (int z = 0; z < dims[2]; z++)
	{
		for (int y = 0; y < dims[1]; y++)
		{
			int offset = z * dims[0] * dims[1] + y * dims[0] + x;
			double val = values->GetValue(offset);
			if (val > 0)
				positiveVal = true;
			if (val < 0)
				negativeVal = true;
			if (positiveVal && negativeVal)
				break;
		}
		if (positiveVal && negativeVal)
			break;
	}

	bool IsoSurf = positiveVal && negativeVal;

	return IsoSurf;
}

bool CMRFSurfaceReconstruction::CheckZSliceIsoSurface(vtkImageData *vol, int z)
{
	int dims[3] = { 0,0,0 };
	vol->GetDimensions(dims);

	vtkDoubleArray *values =
		vtkDoubleArray::SafeDownCast(vol->GetPointData()->GetScalars());

	bool negativeVal = false;
	bool positiveVal = false;
	int zoff = z * dims[0] * dims[1];
	for (int x = 0; x < dims[0]; x++)
	{
		for (int y = 0; y < dims[1]; y++)
		{
			int offset = zoff + y * dims[0] + x;
			double val = values->GetValue(offset);
			if (val > 0)
				positiveVal = true;
			if (val < 0)
				negativeVal = true;
			if (positiveVal && negativeVal)
				break;
		}
		if (positiveVal && negativeVal)
			break;
	}

	bool IsoSurf = positiveVal && negativeVal;

	return IsoSurf;
}


//            0     1     2    3    4   5
// Bounds: (minx, maxx, miny,maxy,minz,maxz)
void CMRFSurfaceReconstruction::GetIsoSurfaceBounds(vtkImageData *vol, double *bounds)
{
	int dims[3] = { 0,0,0 };
	vol->GetDimensions(dims);
	
	double vbounds[6];
	vol->GetBounds(vbounds);
	double spacing[3];
	vol->GetSpacing(spacing);

	// Check x-slices
	bounds[0] = -1;
	bool minXIS = CheckXSliceIsoSurface(vol, 0);
	if (!minXIS) // No iso-surface cross edge
	{
		for (int x = 1; x < dims[0]; x++)
		{
			if (CheckXSliceIsoSurface(vol, x))
			{
				bounds[0] = vbounds[0] + spacing[0] * x;
				break;
			}
		}
	}
	bounds[1] = -1;
	bool maxXIS = CheckXSliceIsoSurface(vol, dims[0] - 1);
	if (!maxXIS) // No iso-surface cross edge
	{
		for (int x = dims[0] - 1; x >= 0; x--)
		{
			if (CheckXSliceIsoSurface(vol, x))
			{
				bounds[1] = vbounds[0] + spacing[0] * x;
				break;
			}
		}
	}

	// Check y-slices
	bounds[2] = -1;
	bool minYIS = CheckYSliceIsoSurface(vol, 0);
	if (!minYIS) // No iso-surface cross edge
	{
		for (int y = 1; y < dims[1]; y++)
		{
			if (CheckYSliceIsoSurface(vol, y))
			{
				bounds[2] = vbounds[2] + spacing[1] * y;
				break;
			}
		}
	}

	bounds[3] = -1;
	bool maxYIS = CheckYSliceIsoSurface(vol, dims[1]-1);
	if (!maxYIS) // No iso-surface cross edge
	{
		for (int y = dims[1] - 1; y >= 0; y--)
		{
			if (CheckYSliceIsoSurface(vol, y))
			{
				bounds[3] = vbounds[2] + spacing[1] * y;
				break;
			}
		}
	}

	// Check z-slices
	bounds[4] = -1;
	bool minZIS = CheckZSliceIsoSurface(vol, 0);
	if (!minZIS) // No iso-surface cross edge
	{
		for (int z = 1; z < dims[2]; z++)
		{
			if (CheckZSliceIsoSurface(vol, z))
			{
				bounds[4] = vbounds[4] + spacing[2] * z;
				break;
			}
		}
	}

	bounds[5] = -1;
	bool maxZIS = CheckZSliceIsoSurface(vol, dims[2] - 1);
	if (!maxZIS) // No iso-surface cross edge
	{
		for (int z = dims[2] - 1; z >= 0; z--)
		{
			if (CheckZSliceIsoSurface(vol, z))
			{
				bounds[5] = vbounds[4] + spacing[2] * z;
				break;
			}
		}
	}

}

int CMRFSurfaceReconstruction::EstimateNewVolumeDimensions( vtkImageData *vol, vtkPolyData *points, int *extent, double &sampleSpacing, double *origin, int maxSize, double MinCellSize )
{
	if (!m_Parms.ShrinkBoundingVolume)
	{
		std::cout << "Estimating new volume without shrinking" << std::endl;
		int size = EstimateNewVolumeDimensionsWithoutShrink(vol, points, extent, sampleSpacing, origin, maxSize, MinCellSize);
		return size;
	}

	std::cout << "Estimating new volume with shrinking" << std::endl;

	// Get distance between voxels in org volume (cellsize)
	double spacing[3];
	vol->GetSpacing(spacing);
	double RealSS = spacing[0];

	// Get dimensions of org volume
//	double VolBounds[6];
//	vol->GetBounds(VolBounds);

	// Get dimensions of point cloud
	double PBounds[6];
	points->GetBounds(PBounds);

	// Get dimensions of possible iso-surface
	double ISbounds[6];
	GetIsoSurfaceBounds(vol, ISbounds);

	// Modify bounds based on previous isosurface
	if (ISbounds[0] != -1 && ISbounds[0] < PBounds[0])  // xmin
		PBounds[0] = ISbounds[0];
	if (ISbounds[1] != -1 && ISbounds[1] > PBounds[1])  // xmax
		PBounds[1] = ISbounds[1];
	if (ISbounds[2] != -1 && ISbounds[2] < PBounds[2])  // Ymin
		PBounds[2] = ISbounds[2];
	if (ISbounds[3] != -1 && ISbounds[3] > PBounds[3])  // Ymax
		PBounds[3] = ISbounds[3];
	if (ISbounds[4] != -1 && ISbounds[4] < PBounds[4])  // Zmin
		PBounds[4] = ISbounds[4];
	if (ISbounds[5] != -1 && ISbounds[5] > PBounds[5])  // Zmax
		PBounds[5] = ISbounds[5];

	// The real size (probably in mm ..) of the point cloud (or the existing iso-surface
	double Xs = PBounds[1] - PBounds[0];
	double Ys = PBounds[3] - PBounds[2];
	double Zs = PBounds[5] - PBounds[4];

	// First try with half the samplespacing
	double ss = RealSS / 2;

	// How many voxels on each side
	int PadSize = m_Parms.PadVoxels;

	int size = maxSize+1;

	// sx, sy, sz are the number of voxels in each dimension
	int sx = 0;
	int sy = 0;
	int sz = 0;

	// ss is the cell side length
	if (ss < MinCellSize)
	{
		ss = MinCellSize;
		std::cout << "cell size hit MinCellSize " << MinCellSize << std::endl;
	}
	do 
	{
		sx = (int)(Xs/ss + 2 * PadSize) + 1;
		sy = (int)(Ys/ss + 2 * PadSize) + 1;
		sz = (int)(Zs/ss + 2 * PadSize) + 1;

		size = sx * sy * sz;

		if (size > maxSize)
		{
			ss *= 1.1;
		}

	} while (size > maxSize);

	sampleSpacing = ss;
	double Xmin = PBounds[0] - PadSize * sampleSpacing;
	double Ymin = PBounds[2] - PadSize * sampleSpacing;
	double Zmin = PBounds[4] - PadSize * sampleSpacing;

	origin[0] = Xmin;
	origin[1] = Ymin;
	origin[2] = Zmin;

	extent[0] = 0;
	extent[1] = sx;
	extent[2] = 0;
	extent[3] = sy;
	extent[4] = 0;
	extent[5] = sz;

//	CheckVolumeIfIsoSurfacePresentAtEdges(vol, points, extent, sampleSpacing, origin, maxSize, MinCellSize);

	return size;
}

bool CMRFSurfaceReconstruction::MultiLevelMRFRegularisation()
{
	std::string iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals.vtk";
	std::string oname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
	std::string oname_mhd			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.mhd";

	std::string onameDistBounds		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.txt";
	std::string onameDistBoundsV	= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Bounds.vtk";
	std::string onameDistDims		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Dimensions.txt";
	std::string onameDistCS			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_CellSize.txt";
	std::string onameVoxel			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_VoxelSize.vtk";
	std::string onameTimings		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_DistanceField_Timings.txt";

	if (MRFDistanceField)
		return true;

	if (m_Parms.UseLeavePatchOut)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_cut.vtk";
	}
	else if (m_Parms.AddNoise)
	{
		iname				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_Orientednormals_Noise.vtk";
	}

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return true;
	}

	if (!NormalData)
	{
		NormalData = vtkExtMisc::MultiReadSurface(iname);
		if (!NormalData)
		{
			std::cerr << "Could not read " << iname << std::endl;
			return false;
		}
	}

	CMultiTimer::MultiTimer().Start("MRF Regularisation", 1);

	std::cout << "Computing initial distance map" << std::endl;

	int minSideVoxels = 8;

	vtkTimerLog *timer = vtkTimerLog::New();
	timer->StartTimer();

	CMultiTimer::MultiTimer().Start("Initial distance map", 1);

	int ComputeMode = VTK_ORIENTEDDISTANCE_FULLVOLUME;
	if (m_Parms.SignedDistance == false)
	{
		ComputeMode = VTK_UNSIGNEDDISTANCE_FULLVOLUME;
	}

	vtkOrientedPointSetDistanceFilter2 *initSD = vtkOrientedPointSetDistanceFilter2::New();
	initSD->SetInputData(NormalData);
	initSD->SetDistanceMode(m_Parms.DistanceMode);
	initSD->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
	initSD->SetCreateWeightVolume(true);
	initSD->SetMinSideVoxels(minSideVoxels);
	initSD->SetNumberOfDistances(m_Parms.NumberOfDistances);
	initSD->SetSearchMode(VTK_SEARCH_MODE_NNEIGHBOURS);
	initSD->SetDistanceMode(m_Parms.DistanceMode);
	initSD->SetComputeMode(ComputeMode);
//	initSD->SetPadVoxels(m_Parms.PadVoxels);
	initSD->SetPadVoxels(5);
	initSD->Update();

	CMultiTimer::MultiTimer().End("Initial distance map");

	if (m_Parms.WriteLevel > 2)
	{
		WriteBoundsFile(initSD->GetOutput(), 0);
		WriteVoxelSizeFile(initSD->GetOutput(), 0);

		std::string init_df_name = m_Parms.OutPutDir + m_Parms.InputNameShort + "_initial_distance_field.mhd";

		std::cout << "Writing initial distance field (MHD)" << std::endl;
		vtkMetaImageWriter* MHDWriter = vtkMetaImageWriter::New();
		MHDWriter->SetInputData(initSD->GetOutput());
		MHDWriter->SetFileName(init_df_name.c_str());
		MHDWriter->Write();
		MHDWriter->Delete();
	}

	vtkImageData *WeigthVolume = initSD->GetWeightVolume();

	std::cout << "Initial MRF regularisation" << std::endl;

	bool BandedICM = false;
	CMultiTimer::MultiTimer().Start("Initial MRF Regularisation", 1);

	vtkDistanceFieldMRFRegularisation *initialMRF = vtkDistanceFieldMRFRegularisation::New();
	initialMRF->SetInputConnection(initSD->GetOutputPort());
	initialMRF->SetWeightVol(WeigthVolume);
	initialMRF->SetIterations(1000);
	initialMRF->SetGlobalBeta(m_Parms.GlobalBeta);
	initialMRF->SetWeightType(m_Parms.WeightMode);
	initialMRF->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
	initialMRF->SetPriorType(m_Parms.PriorType);
	initialMRF->SetOptimisation(m_Parms.Optimisation);
	initialMRF->SetCGTolerance(2);
	initialMRF->SetBandedICM(BandedICM);
	initialMRF->Update();

	CMultiTimer::MultiTimer().End("Initial MRF Regularisation");
	if (m_Parms.WriteLevel > 2)
	{
		WriteDiagnosisFile(initialMRF->GetFinalSplineEnergy(), 0);
	}

	vtkImageData* curVolume = initSD->GetOutput();

	int maxLevels = 20;

	std::vector<vtkImageReslice*> Resamplers(maxLevels, NULL);
	std::vector<vtkOrientedPointSetDistanceFilter2*> SD(maxLevels, NULL);
	std::vector<vtkDistanceFieldMRFRegularisation*> filters(maxLevels, NULL);

//	double minRMS = 0.001;
	double minRMS = 0.000001;
	int iterations = 1000;
	int MaxVolumeSize = m_Parms.MaxVolumeSize;
	double CGTol = 1;
	int id = 0;
	bool stop = false;

	int resDim[3] = {0,0,0};
	curVolume->GetDimensions(resDim);

	double oldSize = resDim[0] * resDim[1] * resDim[2];

	do 
	{
		std::cout << "*****---- Iteration " << id + 1 << " ----*****" << std::endl;
		DWORDLONG freeDoubles = CMemoryUtils::GetAvailableMemoryForDouble();
		std::cout << "System claims space for " << freeDoubles << " doubles" << std::endl;

		int newExtent[6];
		double newOrigin[3];
		double newSampleSpacing;

		int resSize = EstimateNewVolumeDimensions(curVolume, NormalData, newExtent, newSampleSpacing, newOrigin, MaxVolumeSize, m_Parms.MinCellSize);

		if (resSize > MaxVolumeSize || resSize < (oldSize * 1.5) || newSampleSpacing < m_Parms.MinCellSize)
		{
			std::cout << "Limit reached with (" << newExtent[1]  << ", " << newExtent[3] << ", " << newExtent[5] << ") " << resSize;
			if (resSize > MaxVolumeSize)
			{
				cout << " max volume size reached" << std::endl;
			}
			if (resSize < (oldSize * 1.5))
			{
				cout << " increase in volume size to small" << std::endl;
			}
			if (newSampleSpacing < m_Parms.MinCellSize)
			{
				cout << " samplespacing less than minimum spacing" << std::endl;
			}
			stop = true;
		}
		else
		{
			std::cout << "Upsampling field : (" << newExtent[1] << ", " << newExtent[3] << ", " << newExtent[5] << ") " << resSize << " SS " << newSampleSpacing << std::endl;

			Resamplers[id] = vtkImageReslice::New();
			Resamplers[id]->SetInputData(curVolume);
			Resamplers[id]->SetOutputSpacing(newSampleSpacing, newSampleSpacing, newSampleSpacing);
			Resamplers[id]->SetOutputOrigin(newOrigin);
			Resamplers[id]->SetOutputExtent(newExtent);
			Resamplers[id]->SetInterpolationModeToLinear();
			Resamplers[id]->Update();

			oldSize = resSize;
			std::cout << "Created volume: (" << newExtent[1] << ", " << newExtent[3] << ", " << newExtent[5] << ") " << resSize << " SS " << newSampleSpacing << std::endl;

			std::cout << "Calculating banded distance" << std::endl;
			CMultiTimer::MultiTimer().Start("Calculating banded distance", 1);
			SD[id] = vtkOrientedPointSetDistanceFilter2::New();
			SD[id]->SetInputData(NormalData);
			SD[id]->SetInputDistances(Resamplers[id]->GetOutput());
			SD[id]->SetDistanceMode(m_Parms.DistanceMode);
			SD[id]->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
			SD[id]->SetCreateWeightVolume(true);
			SD[id]->SetComputeMode(VTK_ORIENTEDDISTANCE_BAND);
			SD[id]->SetSearchMode(VTK_SEARCH_MODE_NNEIGHBOURS);
			SD[id]->SetDistanceMode(m_Parms.DistanceMode);
			SD[id]->SetNumberOfDistances(m_Parms.NumberOfDistances);
			SD[id]->SetPadVoxels(m_Parms.PadVoxels);
			SD[id]->Update();
			CMultiTimer::MultiTimer().End("Calculating banded distance");

			WeigthVolume = SD[id]->GetWeightVolume();

			if (id > 0 && m_Parms.BandedICM)
			{
				BandedICM = true;
			}

			std::cout << "Regularisation" << std::endl;
			CMultiTimer::MultiTimer().Start("Regularisation", 1);
			filters[id] = vtkDistanceFieldMRFRegularisation::New();
			filters[id]->SetInputConnection(SD[id]->GetOutputPort());
			filters[id]->SetWeightVol(WeigthVolume);
			filters[id]->SetGlobalBeta(m_Parms.GlobalBeta);
			filters[id]->SetWeightType(m_Parms.WeightMode);
			filters[id]->SetWeightValueHigh(m_Parms.LocalWeightMaxDist);
			filters[id]->SetPriorType(m_Parms.PriorType);
			filters[id]->SetOptimisation(m_Parms.Optimisation);
			filters[id]->SetIterations(iterations);
			filters[id]->SetMinRMS(minRMS);
			filters[id]->SetCGTolerance(CGTol);
			filters[id]->SetBandedICM(BandedICM);
			filters[id]->Update();
			CMultiTimer::MultiTimer().End("Regularisation");

			curVolume = filters[id]->GetOutput();

			if (m_Parms.WriteLevel > 2)
			{
				WriteBoundsFile(curVolume, id+1);
				WriteVoxelSizeFile(curVolume, id+1);
				WriteDiagnosisFile(filters[id]->GetFinalSplineEnergy(), id+1);
			}
		}
		if (iterations > 25)
			iterations /= 4;	

		if (CGTol > 0.1)
			CGTol /= 4;

		id++;

		if (id >= maxLevels)
			stop = true;

	} while (!stop);

	timer->StopTimer();
	double ElapsedTime = timer->GetElapsedTime();
	cout << "\nMRF Regularisation took  " << ElapsedTime << " seconds" << std::endl;

	CMultiTimer::MultiTimer().End("MRF Regularisation");

	if (m_Parms.WriteLevel > 2)
	{
		std::ofstream ost4(onameTimings.c_str());
		if (!ost4)
		{
			std::cerr << "Could not write to " << onameTimings << std::endl;
		}
		else
		{
			ost4 << "Multilevel MRF regularisation " << ElapsedTime << " seconds" << std::endl;
		}
	}

	if (m_Parms.WriteLevel > 1)
	{
		std::cout << "Writing MRF distance field" << std::endl;
		vtkXMLImageDataWriter* DSwriter = vtkXMLImageDataWriter::New();
		DSwriter->SetInputData((vtkDataObject*)curVolume);
		DSwriter->SetFileName(oname.c_str());
		DSwriter->Write();
		DSwriter->Delete();

		vtkMetaImageWriter* MHDWriter = vtkMetaImageWriter::New();
		MHDWriter->SetInputData((vtkDataObject*)curVolume);
		MHDWriter->SetFileName(oname_mhd.c_str());
		MHDWriter->Write();
		MHDWriter->Delete();
	}

	if (m_Parms.WriteLevel > 1)
	{
		double bounds[6];
		curVolume->GetBounds(bounds);

		std::ofstream ost(onameDistBounds.c_str());
		if (!ost)
		{
			std::cerr << "Could not write to " << onameDistBounds << std::endl;
		}
		else
		{
			ost << std::setprecision(10);
			for (int i = 0; i < 6; i++)
				ost << bounds[i] << " ";
		}
		double spacing[3];
		curVolume->GetSpacing(spacing);
		double RealSS = spacing[0];


		std::ofstream ost2(onameDistCS.c_str());
		if (!ost2)
		{
			std::cerr << "Could not write to " << onameDistCS << std::endl;
		}
		else
		{
			ost2 << std::setprecision(10);
			ost2 << RealSS << std::endl;
		}
	}

	MRFDistanceField = vtkImageData::New();
	MRFDistanceField->DeepCopy(curVolume);

	initSD->Delete();
	initialMRF->Delete();

	for (int i = 0; i < maxLevels; i++)
	{
		if (Resamplers[i])
			Resamplers[i]->Delete();
		if (SD[i])
			SD[i]->Delete();
		if (filters[i])
			filters[i]->Delete();
	}
	return true;
}

bool CMRFSurfaceReconstruction::OneShotCompute( CParms &parms, vtkPolyData *input, std::string &msg )
{
	m_Parms = parms;

	DeleteAll();

	RawData = vtkPolyData::New();
	RawData->DeepCopy(input);

	if (!InputDataStatistics())
		return false;

	if (!CreateNormals())
		return false;

	if (!MultiLevelMRFRegularisation())
		return false;

	if (!m_Parms.DoComputeSurface)
		return true;

	if (!PolygoniseDistanceMap())
		return false;

	if (!CropSurface())
		return false;

	return true;
}

bool CMRFSurfaceReconstruction::CheckSurfaceConsistency()
{
	std::string iname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface.vtk";
	std::cout << "Checking surface consistency" << std::endl;

	vtkPolyDataReader *surf = vtkExtMisc::SafeReadPolyData(iname);
	if (!surf)
	{
		std::cerr << "Could not read " << iname << std::endl;
		return false;
	}

	if (surf->GetOutput()->GetNumberOfPoints() < 5)
	{
		std::cout << "Iso-surface contains less than 5 points!" << std::endl;
		return false;
	}

	vtkFeatureEdges *features = vtkFeatureEdges::New();
	features->SetInputConnection(surf->GetOutputPort());
	features->SetNonManifoldEdges(1);
	features->SetBoundaryEdges(0);
	features->SetFeatureEdges(0);
	features->SetManifoldEdges(0);
	features->Update();

	int nfeat = features->GetOutput()->GetNumberOfLines();
	features->Delete();

	if (nfeat > 0)
	{
		std::cout << "Surface contains " << nfeat << " non-manifold edges!" << std::endl;
		return false;
	}

	return true;
}


bool CMRFSurfaceReconstruction::Remesh()
{
	std::string inameMRFDist		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFDistanceField.vti";
	std::string inamePD				= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
	std::string outputname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_remeshed.vtk";

	if (m_Parms.InputName != "" && file_exists(outputname.c_str()))
	{
		return true;
	}

	std::cout << "Remeshing surface" << std::endl;

	vtkXMLImageDataReader *signedDist = vtkXMLImageDataReader::New();
	signedDist->SetFileName(inameMRFDist.c_str());
	signedDist->Update();
	if (!signedDist || signedDist->GetOutput()->GetActualMemorySize() == 0)
	{
		std::cerr << "Could not read " << inameMRFDist << std::endl;
		return false;
	}

	vtkPolyDataReader *surf = vtkExtMisc::SafeReadPolyData(inamePD);
	if (!surf)
	{
		std::cerr << "Could not read " << inamePD << std::endl;
		return false;
	}

	vtkImplicitVolume *impvol = vtkImplicitVolume::New();
	impvol->SetVolume(signedDist->GetOutput());

	vtkPolyData *pd = vtkPolyData::New();

	CGELRemeshing remesher;
	remesher.Remesh(impvol, surf->GetOutput(), pd);

	vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
	norms->SplittingOff();
	norms->ConsistencyOn();
	norms->SetInputData(pd);
	norms->Update();

	vtkExtMisc::WritePDVTK(norms->GetOutput(), outputname);

	norms->Delete();
	pd->Delete();
	surf->Delete();
	impvol->Delete();
	signedDist->Delete();

	return true;
}

void CMRFSurfaceReconstruction::CloneTexturesFromVRML()
{
	std::string inameVRML		= m_Parms.InputDir + m_Parms.InputNameShort + ".wrl";
	std::string inameTex		= m_Parms.InputDir + m_Parms.InputNameShort + ".bmp";
	std::string onameTexOrg1    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_orgTextured.vtk";
	std::string onameTexOrg2    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_orgTextured.bmp";

	std::string inamePD			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
	std::string oname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_texture.vtk";
	std::string onameTex		= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_texture.bmp";

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return;
	}
	if (!file_exists(inameVRML.c_str()))
	{
		std::cout << "Could not read " <<inameVRML << std::endl;
		return;
	}

	std::cout << "Cloning texture" << std::endl;

	// Dirty copy
	std::string cmd = "copy " + inameTex + " " + onameTexOrg2;
	CGeneralUtils::FindReplace(cmd, "/", "\\");
	system(cmd.c_str());

	cmd = "copy " + inameTex + " " + onameTex;
	CGeneralUtils::FindReplace(cmd, "/", "\\");
	system(cmd.c_str());

	// First convert VRML to VTK
	vtkVRMLImporter *vrmlin = vtkVRMLImporter::New();
	vrmlin->SetFileName(inameVRML.c_str());
	vrmlin->Update();

	vtkPolyData *pdVRML = (vtkPolyData*)(vrmlin->GetRenderer()->GetActors()->GetLastActor()->GetMapper()->GetInput());

	// Notice that this clean operation collapses all shared vertices
	// This causes the split lines seen between the two textures
	// the alternative is to keep all vertices (num vertices = 3 * num triangles)
	vtkCleanPolyData *clean = vtkCleanPolyData::New();
	clean->SetInputData(pdVRML);

	vtkPolyDataNormals *norm = vtkPolyDataNormals::New();
	norm->SetInputConnection(clean->GetOutputPort());
	norm->ConsistencyOn();
	norm->SplittingOff();
	norm->Update();

	vtkPolyData *refFace = vtkPolyData::New();
	refFace->DeepCopy(norm->GetOutput());
	refFace->GetPointData()->SetScalars(NULL);

	vtkExtMisc::WritePDVTK(refFace, onameTexOrg1);

	vtkPolyDataReader *surf = vtkExtMisc::SafeReadPolyData(inamePD);
	if (!surf)
	{
		std::cerr << "Could not read " << inamePD << std::endl;
		return;
	}

	vtkPolyData *newFace = vtkPolyData::New();
	newFace->DeepCopy(surf->GetOutput());

	vtkDataArray *TexCoords = refFace->GetPointData()->GetTCoords();

	if (!TexCoords)
	{
		std::cerr << "No texture coordinates" << std::endl;
		return;
	}

	vtkPointLocator *ploc = vtkPointLocator::New();
	ploc->SetDataSet(refFace);
	ploc->SetNumberOfPointsPerBucket(1);
	ploc->BuildLocator();

	const int numPts = newFace->GetNumberOfPoints();

	vtkFloatArray *newTCoords = vtkFloatArray::New();
	newTCoords->SetNumberOfComponents(2);
	newTCoords->Allocate(numPts);
	newTCoords->SetName("TCoords");

	double maxDist = 3 * 3;

	double dumU = 0;
	double dumV = 0.996;

	// Do a simple nearest neighbour lookup and copying of (u,v) texture coordinates
	for (int i = 0; i < numPts; i++)
	{
		double p[3];
		double cp[3];

		newFace->GetPoint(i, p);
		int id = ploc->FindClosestPoint(p);
		refFace->GetPoint(id, cp);

		double TC[2];

		if (vtkMath::Distance2BetweenPoints(p, cp) < maxDist)
		{
			double *TCoords = TexCoords->GetTuple(id);
			TC[0] = TCoords[0];
			TC[1] = TCoords[1];
		}
		else
		{
			TC[0] = dumU;
			TC[1] = dumV;
		}

		newTCoords->InsertTuple2(i, TC[0], TC[1]);
	}

	newFace->GetPointData()->SetTCoords(newTCoords);
	newTCoords->Delete();

	vtkExtMisc::WritePDVTK(newFace, oname);

	newFace->Delete();
	refFace->Delete();
	norm->Delete();
	vrmlin->Delete();
	clean->Delete();
	surf->Delete();
}

//
//bool CMRFSurfaceReconstruction::CheckCameraVisibility( C3dMDCamera &cam, double * p,double *n, vtkCellLocator * locator, double &normScore )
//{
//	double *CamPos = cam.GetPosition();
//
////	double distToCam = vtkMath::Distance2BetweenPoints(p, CamPos);
//
//	double v[3];
//	v[0] = CamPos[0] - p[0];
//	v[1] = CamPos[1] - p[1];
//	v[2] = CamPos[2] - p[2];
//
//	normScore = vtkMath::Dot(v, n);
//	if (normScore < 0)
//		return false;
//
//	// push point a little towards the camera ... else the intersection will be at the point
//	double pp[3];
//	pp[0] = p[0] + 0.01 * (CamPos[0] - p[0]);
//	pp[1] = p[1] + 0.01 * (CamPos[1] - p[1]);
//	pp[2] = p[2] + 0.01 * (CamPos[2] - p[2]);
//
//	double tol = 0.0001;
//	double t = 0;
//	double intx[3];
//	double pcoords[3];
//	int subId = 0;
//	int Intersect = locator->IntersectWithLine(CamPos, pp, tol , t, intx, pcoords, subId);
//
//	return (Intersect == 0);
//}

//
//bool CMRFSurfaceReconstruction::CreateCombinedTexture(const std::string &textureName1,const std::string &textureName2, const std::string &TextureOut)
//{
////	vil_image_view<vil_rgb<vxl_byte> > img1 = vil_load(textureName1.c_str());
//	vil_image_view<vxl_byte> img1 = vil_load(textureName1.c_str());
//	if (!img1)
//	{
//		std::cerr << "Could not read " <<  textureName1 << std::endl;
//		return false;
//	}
//	vil_image_view<vxl_byte> img2 = vil_load(textureName2.c_str());
//	if (!img2)
//	{
//		std::cerr << "Could not read " <<  textureName2 << std::endl;
//		return false;
//	}
//
//	if (img1.nplanes() != img2.nplanes() || img1.nj() != img2.nj())
//	{
//		std::cerr << "Texture image dimensions do not fit" << std::endl;
//		return false;
//	}
//
//	int NewX = img1.ni() + img2.ni();
//	int NewY = img1.nj();
//	int nPlanes = img1.nplanes();
//
//	vil_image_view<vxl_byte> imageOut(NewX,NewY,nPlanes);
//
//	// Copy from image 1
//	for (int p = 0; p < nPlanes; ++p)
//	{
//		for (int j = 0; j < NewY; ++j)
//		{
//			for (unsigned int i = 0; i < img1.ni(); ++i)
//			{
//				imageOut(i,j,p)= img1(i,j,p);
//			}
//		}
//	}
//	// Copy from image 2
//	for (int p = 0; p < nPlanes; ++p)
//	{
//		for (int j = 0; j < NewY; ++j)
//		{
//			for (unsigned int i = 0; i < img2.ni(); ++i)
//			{ 
//				imageOut(i + img1.ni(),j,p)= img2(i,j,p);
//			}
//		}
//	}
//
//	// Create area that is used for non-defined triangles
//	for (int p = 0; p < nPlanes; ++p)
//	{
//		for (int j = 0; j < 5; ++j)
//		{
//			for (int i = 0; i < 5; ++i)
//			{
//				imageOut(i, imageOut.nj() -1 - j,p)= 127;
//			}
//		}
//	}
//
//
//	vil_save(imageOut, TextureOut.c_str());
//
//	
//	return true;
//}

//void CMRFSurfaceReconstruction::ReTextureSurface()
//{
//	std::string surfname			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
//	std::string cam1Cname             = m_Parms.CalibrationDir + "calib_1C.tka";
//	std::string cam2Cname             = m_Parms.CalibrationDir + "calib_2C.tka";
////	std::string oname			    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_InCameraSystem.vtk";
//	std::string oname			    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_Retextured.vtk";
//	std::string oname1			    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_WhichCamera.vtk";
//	std::string oname2			    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_WhichCamera_VerticeSplit.vtk";
//	std::string oname3			    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_WhichCamera_Relabeled.vtk";
//	std::string oname4			    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_Retexture_MultiTextures.vtk";
//	std::string TextureName1        = m_Parms.CalibrationDir + "texture_1C.bmp";
//	std::string TextureName2        = m_Parms.CalibrationDir + "texture_2C.bmp";
//	std::string TextureOutName	    = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_Retexture_MultiTextures.png";
//
//
//	//if (m_Parms.InputName != "" && file_exists(oname.c_str()))
//	//{
//	//	return;
//	//}
//
//	if (!MRFSurfaceCropped)
//	{
//		MRFSurfaceCropped = vtkExtMisc::MultiReadSurface(surfname);
//		if (!MRFSurfaceCropped)
//		{
//			std::cerr << "Could not read " << surfname << std::endl;
//			return;
//		}
//	}
//
//	std::cout << "Retexturing" << std::endl;
//	C3dMDCamera cam1C;
//	if (!cam1C.Read(cam1Cname))
//	{
//		std::cerr << "Could not read camera" << std::endl;
//		return;
//	}
//	C3dMDCamera cam2C;
//	if (!cam2C.Read(cam2Cname))
//	{
//		std::cerr << "Could not read camera" << std::endl;
//		return;
//	}
//
//	if (m_Parms.WriteLevel > 0)
//	{
//		if (!CreateCombinedTexture(TextureName1,TextureName2,TextureOutName))
//			return;
//
//	}
//	//vtkPolyData *ShapeInCamSystem = cam1C.TransformShapeIntoCameraSystem(MRFSurfaceCropped);
//
//	//vtkExtMisc::WritePDVTK(ShapeInCamSystem, oname);
//
//	//ShapeInCamSystem->Delete();
//
//	if (m_Parms.WriteLevel > 10)
//	{
//		vtkPolyData *ReTextured = vtkPolyData::New();
//		ReTextured->DeepCopy(MRFSurfaceCropped);
//		ReTextured->GetPointData()->SetScalars(NULL);
//			
//		cam1C.ComputeAllTextureCoordinates(ReTextured);
//
//		vtkExtMisc::WritePDVTK(ReTextured, oname);
//
//		ReTextured->Delete();
//	}
//
//	// Determine which camera that gave which coordinates
//	vtkPolyData *WhichCameraPD = vtkPolyData::New();
//	WhichCameraPD->DeepCopy(MRFSurfaceCropped);
//
//	vtkDoubleArray *scalars = vtkDoubleArray::New();
//	scalars->SetNumberOfComponents(1);
//
//	vtkCellLocator *locator = vtkCellLocator::New();
//	locator->SetDataSet(WhichCameraPD);
//	locator->SetNumberOfCellsPerBucket(1);
//	locator->BuildLocator();
//
//
//	vtkDataArray *normals = WhichCameraPD->GetPointData()->GetNormals();
//	if (!normals)
//	{
//		std::cerr << "No normals" << std::endl;
//		return;
//	}
//
//	scalars->SetNumberOfValues(WhichCameraPD->GetNumberOfPoints());
//	for (int i = 0; i < WhichCameraPD->GetNumberOfPoints(); i++)
//	{
//		scalars->SetValue(i, 0);
//		double p[3];
//		WhichCameraPD->GetPoint(i, p);
//		double n[3];
//		normals->GetTuple(i, n);
//
//		double normScore1C = 0;
//		double normScore2C = 0;
//		bool ViewCam1C = CheckCameraVisibility(cam1C, p, n, locator, normScore1C);
//		bool ViewCam2C = CheckCameraVisibility(cam2C, p, n, locator, normScore2C);
//
//		double ScalVal = 0;
//		if (ViewCam1C && ViewCam2C) 
//		{
//			if (normScore2C > normScore1C)
//			{
//				ScalVal = 4;
//			}
//			else
//			{
//				ScalVal = 3;
//			}
//		}
//		else if (ViewCam1C)
//		{
//			ScalVal = 1;
//		}
//		else if (ViewCam2C)
//		{
//			ScalVal = 2;
//		}
//		scalars->SetValue(i, ScalVal);
//	}
//
//	WhichCameraPD->GetPointData()->SetScalars(scalars);
//	scalars->Delete();
//	locator->Delete();
//
//	if (m_Parms.WriteLevel > 2)
//	{
//		vtkExtMisc::WritePDVTK(WhichCameraPD, oname1);
//	}
//
//	vtkPolyData *splits = SplitVertices(WhichCameraPD);
//
//	if (m_Parms.WriteLevel > 2)
//	{
//		vtkExtMisc::WritePDVTK(splits, oname2);
//	}
//
//	RelabelVisibilityMesh(splits);
//	
//	if (m_Parms.WriteLevel > 2)
//	{
//		vtkExtMisc::WritePDVTK(splits, oname3);
//	}
//
//	vtkPolyData *ReTextured2 = vtkPolyData::New();
//	ReTextured2->DeepCopy(splits);
//
//	cam1C.ComputeAllTextureCoordinatesForThisCamera(ReTextured2, 1);
//	cam2C.ComputeAllTextureCoordinatesForThisCamera(ReTextured2, 2);
//	ReTextured2->GetPointData()->SetScalars(NULL);
//
//	if (m_Parms.WriteLevel > 0)
//	{
//		vtkExtMisc::WritePDVTK(ReTextured2, oname4);
//	}
//
//	ReTextured2->Delete();
//
//
//	WhichCameraPD->Delete();
//	splits->Delete();
//
//}

//void CMRFSurfaceReconstruction::RelabelVisibilityMesh(vtkPolyData *pd)
//{
//	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetScalars());
//	if (!scalars)
//	{
//		std::cerr << "Something wrong with scalars" << std::endl;
//		return;
//	}
//
//	pd->GetPolys()->InitTraversal();
//
//	for (int i = 0; i < pd->GetNumberOfCells(); i++) 
//	{
//		vtkIdType n_pts = -1;
//		const vtkIdType *pts = NULL;
//		pd->GetPolys()->GetNextCell(n_pts, pts); 
//
//		bool OneFound = false;
//		bool TwoFound = false;
//		int Threes = 0;
//		int Fours = 0;
//
//		// copy the corresponding points, scalars and normals across
//		for (int j = 0;j < n_pts; j++) 
//		{
//			double scal = scalars->GetValue(pts[j]);
//
//			if (scal == 1)
//				OneFound = true;
//			else if (scal == 2)
//				TwoFound = true;
//			else if (scal == 3)
//				Threes++;
//			else if (scal == 4)
//				Fours++;
//			//else
//			//	std::cout << "Val: " << scal << std::endl;
//		}
//		double finalVal = 0;
//		if (OneFound && TwoFound)
//			finalVal = 0;
//		else if (OneFound)
//			finalVal = 1;
//		else if (TwoFound)
//			finalVal = 2;
//		else if (Threes > Fours)
//			finalVal = 1;
//		else if (Fours > Threes)
//			finalVal = 2;
//
//		for (int j = 0;j < n_pts; j++) 
//		{
//			scalars->SetValue(pts[j], finalVal);
//		}
//	}
//}


vtkPolyData * CMRFSurfaceReconstruction::SplitVertices(vtkPolyData *pd)
{
	vtkPolyData *newPD = vtkPolyData::New();

	vtkPoints *new_points = vtkPoints::New();

	vtkDataArray *normals = pd->GetPointData()->GetNormals();
	vtkDoubleArray *new_normals = vtkDoubleArray::New();
	new_normals->SetNumberOfComponents(3);
	vtkCellArray *new_polys = vtkCellArray::New();

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetScalars());
	vtkDoubleArray *newscalars = vtkDoubleArray::New();
	newscalars->SetNumberOfComponents(1);

	pd->GetPolys()->InitTraversal();

	for (int i = 0; i < pd->GetNumberOfCells(); i++) 
	{
		vtkIdType n_pts = -1;
		const vtkIdType *pts = NULL;
		pd->GetPolys()->GetNextCell(n_pts, pts); 

		// Done 26-11-2021: Might be buggy
		vtkIdType* new_pts = new vtkIdType[n_pts];

		// copy the corresponding points, scalars and normals across
		for (int j = 0;j < n_pts; j++) 
		{
			
			if (normals)
			{
				new_normals->InsertNextTuple(normals->GetTuple(pts[j]));
			}
			if (scalars)
			{
				newscalars->InsertNextValue(scalars->GetValue(pts[j]));
			}
			// copy the vertex into the new structure and update
			// the vertex index in the polys structure (pts is a pointer into it)
			new_pts[j] = new_points->InsertNextPoint(pd->GetPoints()->GetPoint(pts[j]));
			
		}
		// copy this poly (pointing at the new points) into the new polys list 
		new_polys->InsertNextCell(n_pts,new_pts);

		delete[] new_pts;
	}

	// use the new structures for the output
	newPD->SetPoints(new_points);
	newPD->SetPolys(new_polys);
	if (normals)
	{
		newPD->GetPointData()->SetNormals(new_normals);
	}
	if (scalars)
	{
		newPD->GetPointData()->SetScalars(newscalars);
	}
	newPD->Squeeze();

	new_points->Delete();
	new_polys->Delete();
	new_normals->Delete();
	newscalars->Delete();
	
	return newPD;
}


bool CMRFSurfaceReconstruction::OneShotComputePCANormals(CParms &parms, vtkPolyData *input, std::string &msg )
{
	m_Parms = parms;

	DeleteAll();

	RawData = vtkPolyData::New();
	RawData->DeepCopy(input);

	if (m_Parms.AdaptiveParameters)
	{
		if (!InputDataStatistics())
		{
			msg = "Could not calculate input statistics";
			return false;
		}
	}

	if (!CreateNormals())
	{
		msg = "Could not calculate normals";
		return false;
	}
	return true;
}

bool CMRFSurfaceReconstruction::OneShotComputeDistanceField( CParms &parms, vtkPolyData *input, std::string &msg )
{
	m_Parms = parms;

	DeleteAll();

	RawData = vtkPolyData::New();
	RawData->DeepCopy(input);

	if (!InputDataStatistics())
		return false;

	if (!CreateNormals())
		return false;

	if (!MultiLevelMRFRegularisation())
		return false;

	return true;
}

int CMRFSurfaceReconstruction::PolygoniseAndRemesh( vtkImageData * SDF, vtkImplicitVolume * impvol, vtkPolyData * output, int Polygoniser, double CellSizeFactor, 
	int remesh, vtkPolyData *targetEdgeLengthPD, bool FixNonManifold, double IsoValue, std::string &msg )
{
 	//CGELRemeshing rem;
 	//rem.Polygonise(SDF, impvol);

	double space[3];
	SDF->GetSpacing(space);

	double RealSS = space[0];

	vtkPolyData *pd = vtkPolyData::New();
	if (Polygoniser == 0)
	{
		CMultiTimer::MultiTimer().Start("Bloomenthal", 1);
		double cellSize = RealSS / CellSizeFactor;

		int resDim[3];
		SDF->GetDimensions(resDim);
		int MaxCells = std::max(resDim[0], std::max(resDim[1], resDim[2])) * 2;

		bool useTetra = false;
		std::cout << "\nBloomenthal with cell size " << cellSize << std::endl;

		double bounds[6];
		SDF->GetBounds(bounds);

		double start[3];
		start[0] = (bounds[1] + bounds[0]) / 2; 
		start[1] = (bounds[3] + bounds[2]) / 2;
		start[2] = (bounds[5] + bounds[4]) / 2;

		try
		{
			bounds[0] += 1.1 * cellSize;
			bounds[1] -= 1.1 * cellSize;
			bounds[2] += 1.1 * cellSize;
			bounds[3] -= 1.1 * cellSize;
			bounds[4] += 1.1 * cellSize;
			bounds[5] -= 1.1 * cellSize;

			Geometry::BloomenthalPolygonizer bloomenthal(impvol, cellSize, MaxCells, bounds, useTetra);
			bloomenthal.march(start[0], start[1], start[2]);
			bloomenthal.ExportToPolyData(pd);
		}
		catch (std::string ermsg)
		{
			std::cerr << "Bloomenthal fail " << ermsg << std::endl;
			msg = "Bloomenthal fail " + ermsg;
			return 0;
		}
		CMultiTimer::MultiTimer().End("Bloomenthal");
	}
	else
	{
		CMultiTimer::MultiTimer().Start("Marching cubes", 1);

		std::cout << "Marching cubes" << std::endl;
		int resDim[3];
		SDF->GetDimensions(resDim);

		// std::cout << "DEBUG: cellsizefactor " << CellSizeFactor << std::endl;
		// Removed cellsizefactor to avoid the sample function creating weird artifacts at the edges.
		int sampDims[3];
		//sampDims[0] = int(resDim[0] * CellSizeFactor);
		//sampDims[1] = int(resDim[1] * CellSizeFactor);
		//sampDims[2] = int(resDim[2] * CellSizeFactor);

		sampDims[0] = int(resDim[0]);
		sampDims[1] = int(resDim[1]);
		sampDims[2] = int(resDim[2]);

		double bounds[6];
		SDF->GetBounds(bounds);

//		std::cout << "DEBUG: before sample" << std::endl;
		vtkSampleFunction *sample = vtkSampleFunction::New();
		sample->SetSampleDimensions(sampDims);
		sample->SetModelBounds(SDF->GetBounds());
		sample->SetImplicitFunction(impvol);
		sample->ComputeNormalsOff();
		sample->SetCapping(false);
//		sample->DebugOn();
		// sample->PrintSelf(std::cout,vtkIndent());
		sample->Update();
//		std::cout << "DEBUG: after sample" << std::endl;


		bool debug = false; 
		if (debug)
		{
			std::cout << "Writing distance field (MHD)" << std::endl;
			std::string dfName = "D:\\Data\\test\\debugout\\Upsampled_DistanceField.mhd";
			vtkMetaImageWriter* MHDWriter = vtkMetaImageWriter::New();
			MHDWriter->SetInputData(sample->GetOutput());
			MHDWriter->SetFileName(dfName.c_str());
			MHDWriter->Write();
			MHDWriter->Delete();
		}


		vtkContourFilter *contour = vtkContourFilter::New();
		contour->SetInputConnection(sample->GetOutputPort());
		contour->GenerateValues(1, IsoValue, IsoValue);
		contour->Update();
//		std::cout << "DEBUG: after contour" << std::endl;

		pd->DeepCopy(contour->GetOutput());

//		std::cout << "DEBUG: after deepcopy" << std::endl;

		contour->Delete();
		sample->Delete();
		CMultiTimer::MultiTimer().End("Marching cubes");
	}

	if (remesh == 0)
	{
		output->DeepCopy(pd);
		pd->Delete();
		return 1;
	}

	if (remesh == 1)
	{
		CMultiTimer::MultiTimer().Start("CheckSurfaceConsistency", 1);
		if (!vtkExtMisc::CheckSurfaceConsistency(pd, msg))
		{
			CMultiTimer::MultiTimer().End("CheckSurfaceConsistency");
			std::cerr << msg << std::endl;
			if (pd->GetNumberOfPoints() > 3 && FixNonManifold)
			{
				std::cout << "Trying to fix non-manifold edges" << std::endl;

				vtkRemoveNonManifoldCells* remover = vtkRemoveNonManifoldCells::New();
				remover->SetInputData(pd);
				remover->Update();

				if (!vtkExtMisc::CheckSurfaceConsistency(remover->GetOutput(), msg))
				{
					std::cout << "Could not remove non-manifold edges" << std::endl;
					output->DeepCopy(pd);
					pd->Delete();
					remover->Delete();
					return 2;
				}
				std::cout << "Succesfully removed non-manifold edges. Now remeshing" << std::endl;

				vtkPolyData* pd2 = vtkPolyData::New();

				CGELRemeshing remesher;
				bool RemeshRes = true;
				if (targetEdgeLengthPD != NULL)
				{
					// Lets just check if any target edge lenghts exists
					vtkDoubleArray* scalars = vtkDoubleArray::SafeDownCast(targetEdgeLengthPD->GetPointData()->GetScalars());
					if (scalars)
					{
						RemeshRes = remesher.RemeshWithCurvature(impvol, remover->GetOutput(), targetEdgeLengthPD, pd2);
					}
					else
					{
						std::cout << "No target edge length scalars exist - normal uniform remeshing" << std::endl;
						RemeshRes = remesher.Remesh(impvol, remover->GetOutput(), pd2);
					}

				}
				else
				{
					RemeshRes = remesher.Remesh(impvol, remover->GetOutput(), pd2);
				}
				if (!RemeshRes)
				{
					pd2->Delete();
					std::cout << "Raw iso-surface is returned since GEL can not remesh " << std::endl;
					output->DeepCopy(pd);
					pd->Delete();
					remover->Delete();
					return 2;
				}

				CMultiTimer::MultiTimer().Start("Computing normals", 1);
				vtkPolyDataNormals* norms = vtkPolyDataNormals::New();
				norms->SplittingOff();
				norms->ConsistencyOn();
				norms->SetInputData(pd2);
				norms->Update();
				CMultiTimer::MultiTimer().End("Computing normals");

				output->DeepCopy(norms->GetOutput());
				pd2->Delete();
				norms->Delete();
				remover->Delete();
				pd->Delete();
				msg = "Remeshed surface after removal of non-manifold edges";
				return 2;
			}
		}
		else
		{
			CMultiTimer::MultiTimer().End("CheckSurfaceConsistency");
			vtkPolyData* pd2 = vtkPolyData::New();

			CGELRemeshing remesher;
			bool RemeshRes = true;
			if (targetEdgeLengthPD != NULL)
			{
				// Lets just check if any target edge lenghts exists
				vtkDoubleArray* scalars = vtkDoubleArray::SafeDownCast(targetEdgeLengthPD->GetPointData()->GetScalars());
				if (scalars)
				{
					RemeshRes = remesher.RemeshWithCurvature(impvol, pd, targetEdgeLengthPD, pd2);
				}
				else
				{
					std::cout << "No target edge length scalars exist - normal uniform remeshing" << std::endl;
					RemeshRes = remesher.Remesh(impvol, pd, pd2);
				}
			}
			else
			{
				double edgeFactor = 1;
				if (CellSizeFactor)
					edgeFactor = edgeFactor / CellSizeFactor;
				RemeshRes = remesher.Remesh(impvol, pd, pd2, -1, edgeFactor);
			}
			if (!RemeshRes)
			{
				pd2->Delete();
				std::cout << "Raw iso-surface is returned since GEL can not remesh " << std::endl;
				output->DeepCopy(pd);
				pd->Delete();
				return 2;
			}

			vtkPolyDataNormals* norms = vtkPolyDataNormals::New();
			norms->SplittingOff();
			norms->SetInputData(pd2);
			norms->ConsistencyOn();
			norms->Update();

			output->DeepCopy(norms->GetOutput());
			pd2->Delete();
			norms->Delete();
		}
		pd->Delete();
		return 1;
	}
	if (remesh == 2)
	{
		std::cout << "Remeshing unsigned distance field using bisection based projection" << std::endl;

		vtkPolyData* pd2 = vtkPolyData::New();

		CGELRemeshing remesher;
		// bool RemeshRes = remesher.ProjectOnly(impvol, pd, pd2);
		bool RemeshRes = remesher.RemeshWithBisectionProjection(impvol, pd, pd2);
		if (!RemeshRes)
		{
			pd2->Delete();
			std::cout << "Raw iso-surface is returned since GEL can not remesh " << std::endl;
			output->DeepCopy(pd);
			pd->Delete();
			return 2;
		}
		else
		{
			output->DeepCopy(pd2);
			pd2->Delete();
			pd->Delete();
			return 1;
		}
	}
	if (remesh == 3)
	{
		std::cout << "Remeshing unsigned distance field using point based projection" << std::endl;

		vtkPolyData* pd2 = vtkPolyData::New();

		CGELRemeshing remesher;
		// bool RemeshRes = remesher.ProjectOnly(impvol, pd, pd2);
		bool RemeshRes = remesher.RemeshWithMultiPointProjection(impvol, pd, pd2);
		if (!RemeshRes)
		{
			pd2->Delete();
			std::cout << "Raw iso-surface is returned since GEL can not remesh " << std::endl;
			output->DeepCopy(pd);
			pd->Delete();
			return 2;
		}
		else
		{
			output->DeepCopy(pd2);
			pd2->Delete();
			pd->Delete();
			return 1;
		}
	}
	return 1;
}

bool CMRFSurfaceReconstruction::SurfaceStatistics()
{
	std::string inamePD			= m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped.vtk";
	std::string oname           = m_Parms.OutPutDir + m_Parms.InputNameShort + "_MRFSurface_cropped_statistics.txt";

	if (m_Parms.InputName != "" && file_exists(oname.c_str()))
	{
		return true;
	}

	std::cout << "Calculating edge statistics" << std::endl;

	vtkPolyDataReader *reader = vtkExtMisc::SafeReadPolyData(inamePD);
	if (!reader)
	{
		std::cerr << "Could not read " << inamePD << std::endl;
		return false;
	}

	std::string msg;
	if (!vtkExtMisc::CheckSurfaceConsistency(reader->GetOutput(), msg))
	{
		std::cerr << msg << std::endl;
		reader->Delete();
		return false;
	}

	std::vector<double> lengths;
	CGELRemeshing::ComputeEdgeStatistics(reader->GetOutput(), lengths);


	int NEdges = lengths.size();
	double frac05 = 0;
	double frac95 = 0;
	double frac99 = 0;
	double median = 0;
	double mean = 0;
	double sdev = 0;
	double minVal = 0;
	double maxVal = 0;
	CGeneralUtils::MeanAndSdev(lengths, mean, sdev);
	CGeneralUtils::MinMax(lengths, minVal, maxVal);
	CGeneralUtils::Median(lengths, 0.05, frac05);
	CGeneralUtils::Median(lengths, 0.95, frac95);
	CGeneralUtils::Median(lengths, 0.99, frac99);
	CGeneralUtils::Median(lengths, 0.50, median);
	double RMS = CGeneralUtils::RMS(lengths);

	std::ofstream ost(oname.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << oname << std::endl;
	}
	else
	{
		ost << "Statistics of mesh edge lengths" << std::endl;
		ost << "Number of edges " << NEdges << std::endl;
		ost << " Min: "  << minVal << std::endl;
		ost << " Max: "  << maxVal << std::endl;
		ost << " Average: "  << mean << std::endl;
		ost << " Sdev: "  << sdev << std::endl;
		ost << " Median: "  << median << std::endl;
		ost << " RMS: "  << RMS << std::endl;
		ost << " 5% fractile: "  << frac05 << std::endl;
		ost << " 95% fractile: "  << frac95 << std::endl;
		ost << " 99% fractile: "  << frac99 << std::endl;
	}
	reader->Delete();

	return true;
}

vtkPolyData * CMRFSurfaceReconstruction::GetNormalData() const
{
	return NormalData;
}

vtkPolyData * CMRFSurfaceReconstruction::GetMRFSurface() const
{
	return MRFSurface;
}

vtkPolyData * CMRFSurfaceReconstruction::GetMRFSurfaceCropped() const
{
	return MRFSurfaceCropped;
}

vtkImageData * CMRFSurfaceReconstruction::GetMRFDistanceField() const
{
	return MRFDistanceField;
}

CMRFSurfaceReconstruction::CParms CMRFSurfaceReconstruction::GetUpdatedParms() const
{
	return m_Parms;
}

vtkPolyData * CMRFSurfaceReconstruction::GetMRFSurfaceAgressivelyCropped() const
{
	return MRFSurfaceAggresivelyCropped;
}

CMRFSurfaceReconstruction::CParms::CParms()
{
	InputType = 1;
	EstimateNormals = true;
	SignedDistance = true;
	ConnectedComponentAnalysis = true;
	UseReferenceNormalsInVoting = false;
	WriteLevel = 5;
	ScaleInput = false;
	MaxPlaneDistance = 1.0;
	NormalRadius = 2.0;
	NormalRadiusFactor = 2.5;
	NormalSearchMode = 0;
	NormalPointsPerNormal = 5;
	Max3rdEigenval = 10.0;
	MinCliqueSize = 500;
	MinCliquePercent = 0.01;
	CliqueNeighbourDistance = 2.5;
	MinCellSize = 2.5;
	LargestCliqueOnly = false;
//	RotateToMinVolume = false;
	SampleSpace = -1;
	SampleFactor = 2;
	kappa = 3;
	LocalWeightMaxDist = 3;
//	EnlargePercent = 50;
	PadVoxels = 5;
	ShrinkBoundingVolume = true;
	DistanceMode = 2;
	NumberOfDistances = 5;
	PriorType = 0;
	UseLocalWeights = false;
	MaxIterations = 1000;
	GlobalBeta = 0.1;
	CellSizeFactor = 2.0;
	BandedICM = true;
	Polygoniser = 1;
	WeightMode = 1;
	//OhtakeOptimise = false;
	UseLeavePatchOut = false;
	AddNoise = false;
	NoiseNature = 0;
	MultiLevel = false;
	MaxVolumeSize = 10000000;
	Remesh = true;
	RemeshToTargetEdgeLengths = false;
	CloneTextures = false;
	Optimisation = 1;
	AdaptiveParameters = true;
	AggresiveCrop = false;
	DoComputeSurface = 1;
	IsoValue = 0;
}
/*
0 - m1, m2 
1 3dMD scanner m1, m3 
2 - (m1), (m3) 
3 3Shape earscan, Stanford models m4, (m1), (m3) 
4 - (m1), m3 
5 Many 3D models (m1), (m3) 

*/

void CMRFSurfaceReconstruction::CParms::SetInputType( int type )
{
	InputType = type;
	if (type == 0)
	{
		EstimateNormals = true;
		ConnectedComponentAnalysis = true;
		UseReferenceNormalsInVoting = false;
	}
	else if (type == 1)
	{
		EstimateNormals = true;
		ConnectedComponentAnalysis = true;
		UseReferenceNormalsInVoting = true;
	}
	else if (type == 2)
	{
		EstimateNormals = false;
		ConnectedComponentAnalysis = false;
		UseReferenceNormalsInVoting = true;
	}
	else if (type == 3)
	{
		EstimateNormals = true;
		ConnectedComponentAnalysis = true;
		UseReferenceNormalsInVoting = true;
	}
	else if (type == 4)
	{
		EstimateNormals = true;
		ConnectedComponentAnalysis = true;
		UseReferenceNormalsInVoting = true;
	}
	else if (type == 5)
	{
		EstimateNormals = false;
		ConnectedComponentAnalysis = false;
		UseReferenceNormalsInVoting = true;
	}
}

void CMRFSurfaceReconstruction::CParms::DumpToFile(const std::string &name)
{
	std::string fname = OutPutDir + InputNameShort + "_parameters.txt";

	if (name != "")
	{
		fname = name;
	}

	std::ofstream fost(fname.c_str());
	if (!fost)
	{
		std::cerr << "Could not write to " << fname << std::endl;
		return;
	}
	std::string doneByAdaptive = "";
	if (AdaptiveParameters)
	{
		doneByAdaptive = " (estimated)";
	}

	fost << "Output dir "  << OutPutDir << std::endl;
	fost << "Input dir "  << InputDir << std::endl;
	fost << "Input name " << InputName << std::endl;
	fost << "Short input name " << InputNameShort << std::endl;
	fost << "WriteLevel " << WriteLevel << std::endl;
	fost << "Calibration Dir " << CalibrationDir << std::endl;
	fost << "InputType " << InputType << ":";
	if (InputType == 0 ) 
		fost << " Points with no normals" << std::endl;
	if (InputType == 1 ) 
		fost << " Points with approximate normals. Normals can be estimated roughly by an acquisition device" << std::endl;
	if (InputType == 2 ) 
		fost << " Points with good normals. Points can be noise/outliers but still have good normal estimates" << std::endl;
	if (InputType == 3 ) 
		fost << " Mesh with connectivity but no normals. STL files for example." << std::endl;
	if (InputType == 4 ) 
		fost << " Mesh with approximate normals. Mesh with connectivity and rough estimates of normal directions." << std::endl;
	if (InputType == 5 ) 
		fost << " Mesh with good normals. Mesh with connectivity and well defined normals" << std::endl;
	fost << "Signed Distance " << SignedDistance << std::endl;
	fost << "EstimateNormals " << EstimateNormals << std::endl;
	fost << "ConnectedComponentAnalysis " << ConnectedComponentAnalysis << std::endl;
	fost << "UseReferenceNormalsInVoting " << UseReferenceNormalsInVoting << std::endl;
	fost << "ScaleInput " << ScaleInput << std::endl;
	fost << "Normal radius " << NormalRadius << doneByAdaptive <<  std::endl;
	fost << "Normal radius factor " << NormalRadiusFactor << std::endl;
	fost << "NormalSearchMode " << NormalSearchMode;
	if (NormalSearchMode == 0 ) 
		fost << ": search in the search radius for points" << std::endl;
	else
		fost << ": search for NormalPointsPerNormal" << std::endl;
	fost << "NormalPointsPerNormal " << NormalPointsPerNormal << std::endl;
	fost << "Max3rdEigenval " << Max3rdEigenval << std::endl;
	fost << "MaxPlaneDistance " << MaxPlaneDistance << doneByAdaptive << std::endl;
	fost << "MinCliquePercent " << MinCliquePercent << std::endl;
	fost << "MinCliqueSize " << MinCliqueSize << doneByAdaptive << std::endl;
	fost << "MinCellSize " << MinCellSize << doneByAdaptive << std::endl;
	fost << "LargestCliqueOnly " << LargestCliqueOnly << std::endl;
	fost << "CliqueNeighbourDistance " << CliqueNeighbourDistance << doneByAdaptive << std::endl;
	fost << "SampleSpace " << SampleSpace<< std::endl;
	fost << "kappa " << kappa << std::endl;
	fost << "LocalWeightMaxDist " << LocalWeightMaxDist << doneByAdaptive << std::endl;
	fost << "SampleFactor " << SampleFactor<< std::endl;
//	fost << "EnlargePercent " << EnlargePercent << std::endl;
	fost << "PadVoxels " << PadVoxels << std::endl;
	fost << "ShrinkBoundingVolume " << ShrinkBoundingVolume << std::endl;
	fost << "BandedICM " << BandedICM << std::endl;
	fost << "DistanceMode " << DistanceMode << ":";
	if (DistanceMode == 0)
		fost << " Nearest point Euclidean distance" << std::endl;
	if (DistanceMode == 1)
		fost << " Nearest point projected distance" << std::endl;
	if (DistanceMode == 2)
		fost << " Median of nearest points projected distance (L1 norm)" << std::endl;
	if (DistanceMode == 3)
		fost << " Average of nearest points projected distance (L2 norm)" << std::endl;
	fost << "Polygoniser " << Polygoniser << std::endl;
	fost << "NumberOfDistances " << NumberOfDistances << std::endl;
	fost << "PriorType " << PriorType << ":";
	if (PriorType == 0)
		fost << " Difference of Laplacians" << std::endl;
	if (PriorType == 1)
		fost << " Difference of voxel values" << std::endl;
	fost << "UseLocalWeights " << UseLocalWeights<< std::endl;
	fost << "WeightMode " << WeightMode<< std::endl;
	fost << "MaxIterations " << MaxIterations<< std::endl;
	fost << "GlobalBeta " << GlobalBeta<< std::endl;
	fost << "CellSizeFactor " << CellSizeFactor<< std::endl;
	fost << "UseLeavePatchOut " << UseLeavePatchOut << std::endl;
	fost << "AddNoise " << AddNoise << std::endl;
	fost << "NoiseNature " << NoiseNature << std::endl;
	fost << "MultiLevel " << MultiLevel << std::endl;
	fost << "MaxVolumeSize " << MaxVolumeSize << std::endl;
	fost << "Remesh " << Remesh << std::endl;
	fost << "Remesh to target edge lengths " << RemeshToTargetEdgeLengths << std::endl;
	fost << "CloneTextures " << CloneTextures << std::endl;
	fost << "Optimisation " << Optimisation << ":";
	if (Optimisation == 0)
		fost << " ICM" << std::endl;
	if (Optimisation == 1)
		fost << " Conjugate Gradient" << std::endl;
	if (Optimisation == 2)
		fost << " Sparse Cholesky" << std::endl;
	fost << "AdaptiveParameters " << AdaptiveParameters << std::endl;
	fost << "AggresiveCrop " << AggresiveCrop << std::endl;
}
