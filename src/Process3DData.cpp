#include "Process3DData.h"


#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkButterflySubdivisionFilter.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricClustering.h>
#include <vtkQuadricDecimation.h>
#include <vtkFeatureEdges.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkConstrainedICPTransform.h>
#include "vtkExtMisc.h"
#include "vtkRemoveUnusedPolyDataPoints.h"

#include "vtkPolyDataDifference.h"

CProcess3DData::CProcess3DData()
{
}

CProcess3DData::~CProcess3DData()
{
}


CICPParameters::CICPParameters()
{
	SourceID = 0;
	TargetID = 1;
	ApplyTransformTo = eSource;
	transSurfaceID = 1;
	SignedDistance = false;
	MatchCentroids = false; 
	MaxLandmarks = 1000;
	MaxIterations = 100;
	MaxMeanDistance = 0.001;
	ExcludeEdgePoints = false;
	CalculateRegErrors = true;
	source = NULL;
	target = NULL;
	sourcePreTransform = NULL;
	targetPreTransform = NULL;
	transformedSurface = NULL;
//	surfaceErrors = NULL;
}

void CProcess3DData::DoICP(CICPParameters& parms)
{
	vtkTransform *sourceTrans = vtkTransform::New();
	 sourceTrans->SetMatrix(parms.sourcePreTransform);

	vtkTransformPolyDataFilter *ptransSource = vtkTransformPolyDataFilter::New();
	 ptransSource->SetInputData(parms.source);
	 ptransSource->SetTransform(sourceTrans);
	 ptransSource->Update();

	 vtkTransform *targetTrans = vtkTransform::New();
	  targetTrans->SetMatrix(parms.targetPreTransform);

	vtkTransformPolyDataFilter *ptransTarget = vtkTransformPolyDataFilter::New();
	 ptransTarget->SetInputData(parms.target);
	 ptransTarget->SetTransform(targetTrans);
	 ptransTarget->Update();


	vtkTransform *moveTrans = vtkTransform::New();
	moveTrans->SetMatrix(parms.moveSurfacePreTransform);

	vtkTransformPolyDataFilter *ptransmove = vtkTransformPolyDataFilter::New();
	ptransmove->SetInputData(parms.move);
	ptransmove->SetTransform(moveTrans);
	ptransmove->Update();

	if (parms.ExcludeEdgePoints)
	{
		vtkConstrainedICPTransform *ICP = vtkConstrainedICPTransform::New();
		ICP->SetSource(ptransSource->GetOutput());
		ICP->SetTarget(ptransTarget->GetOutput());
		ICP->CheckMeanDistanceOn();
		ICP->SetMaximumMeanDistance(parms.MaxMeanDistance);
		ICP->SetMaximumNumberOfIterations(parms.MaxIterations);
		ICP->SetMaximumNumberOfLandmarks(parms.MaxLandmarks);
		ICP->SetStartByMatchingCentroids(parms.MatchCentroids);
		ICP->GetLandmarkTransform()->SetModeToRigidBody();
		ICP->Update();		

		vtkTransformPolyDataFilter *trans = vtkTransformPolyDataFilter::New();
		trans->SetInputConnection(ptransmove->GetOutputPort());
		trans->SetTransform(ICP);
		trans->Update();

		if (parms.CalculateRegErrors)
		{
			vtkPolyDataDifference *PDDiff = vtkPolyDataDifference::New();
			PDDiff->SetSignedDistance(parms.SignedDistance);
			PDDiff->SetExcludeEdgePoints(true);

			PDDiff->SetInputConnection(trans->GetOutputPort());
			PDDiff->SetTargetData(ptransTarget->GetOutput());
			PDDiff->Update();

			parms.transformedSurface->DeepCopy(PDDiff->GetOutput());
			PDDiff->Delete();
		}
		else
		{
			parms.transformedSurface->DeepCopy(trans->GetOutput());
		}
		ICP->Delete();
		trans->Delete();
	}
	else
	{
		vtkIterativeClosestPointTransform *ICP = vtkIterativeClosestPointTransform::New();
		ICP->SetSource(ptransSource->GetOutput());
		ICP->SetTarget(ptransTarget->GetOutput());
		ICP->CheckMeanDistanceOn();
		ICP->SetMaximumMeanDistance(parms.MaxMeanDistance);
		ICP->SetMaximumNumberOfIterations(parms.MaxIterations);
		ICP->SetMaximumNumberOfLandmarks(parms.MaxLandmarks);
		ICP->SetStartByMatchingCentroids(parms.MatchCentroids);
		ICP->GetLandmarkTransform()->SetModeToRigidBody();
		ICP->Update();		

		vtkTransformPolyDataFilter *trans = vtkTransformPolyDataFilter::New();
		trans->SetInputConnection(ptransmove->GetOutputPort());
		trans->SetTransform(ICP);
		trans->Update();

		if (parms.CalculateRegErrors)
		{
			vtkPolyDataDifference *PDDiff = vtkPolyDataDifference::New();
			PDDiff->SetSignedDistance(parms.SignedDistance);

			PDDiff->SetInputConnection(trans->GetOutputPort());
			PDDiff->SetTargetData(ptransTarget->GetOutput());
			PDDiff->Update();

			parms.transformedSurface->DeepCopy(PDDiff->GetOutput());
			PDDiff->Delete();
		}
		else
		{
			parms.transformedSurface->DeepCopy(trans->GetOutput());
		}
		ICP->Delete();
		trans->Delete();
	}
	ptransTarget->Delete();
	targetTrans->Delete();
	ptransSource->Delete();
	sourceTrans->Delete();
	ptransmove->Delete();
	moveTrans->Delete();
}

void CProcess3DData::DoSubdivideSurface(vtkPolyData *source, int subdivisions, vtkPolyData *subdivided, 
			int subdivtype)
{
	vtkTriangleFilter *tri = vtkTriangleFilter::New();
	tri->SetInputData(source);
	tri->Update();
	if (subdivtype == 0)
	{
		vtkLinearSubdivisionFilter *subdiv = vtkLinearSubdivisionFilter::New();
		subdiv->SetInputConnection(tri->GetOutputPort());
		subdiv->SetNumberOfSubdivisions(subdivisions);
		subdiv->Update();
		subdivided->DeepCopy(subdiv->GetOutput());
		subdiv->Delete();
	}
	else if (subdivtype == 1)
	{
		vtkButterflySubdivisionFilter *subdiv = vtkButterflySubdivisionFilter::New();
		subdiv->SetInputConnection(tri->GetOutputPort());
		subdiv->SetNumberOfSubdivisions(subdivisions);
		subdiv->Update();
		subdivided->DeepCopy(subdiv->GetOutput());
		subdiv->Delete();
	}
	else if (subdivtype == 2)
	{
		vtkLoopSubdivisionFilter *subdiv = vtkLoopSubdivisionFilter::New();
		subdiv->SetInputConnection(tri->GetOutputPort());
		subdiv->SetNumberOfSubdivisions(subdivisions);
		subdiv->Update();
		subdivided->DeepCopy(subdiv->GetOutput());
		subdiv->Delete();
	}
	tri->Delete();
}

void CProcess3DData::DoDecimateSurface(vtkPolyData* source, vtkPolyData* decimated,
	int decimtype, float decimfactor, bool preservetopology)
{
	if (decimtype == 2)
	{
		vtkNew<vtkQuadricClustering> decimate;
		decimate->SetInputData(source);
		// decimate->UseFeatureEdgesOn();
		decimate->Update();
		decimated->DeepCopy(decimate->GetOutput());
	}
	else
	{
		vtkTriangleFilter* tri = vtkTriangleFilter::New();
		tri->SetInputData(source);
		tri->Update();
		if (decimtype == 0)
		{
			vtkDecimatePro* decim = vtkDecimatePro::New();
			decim->SetInputConnection(tri->GetOutputPort());
			decim->SetPreserveTopology(preservetopology);
			decim->SetTargetReduction(decimfactor);
			decim->Update();
			decimated->DeepCopy(decim->GetOutput());
			decim->Delete();
		}
		if (decimtype == 1)
		{
			vtkQuadricDecimation* decim = vtkQuadricDecimation::New();
			decim->SetInputConnection(tri->GetOutputPort());
			decim->SetTargetReduction(decimfactor);
			decim->Update();
			decimated->DeepCopy(decim->GetOutput());
			decim->Delete();
		}
		tri->Delete();
	}
}

// #include <vtkConstrainedSmoothingFilter.h>

void CProcess3DData::DoSmoothSurface(vtkPolyData *source, vtkPolyData *smooth, int smoothtype, int NumIt, double RelaxFactor,
			bool BoundarySmooth, bool FeatureEdgeSmooth, double FeatureAngle, bool GenerateErrScal)
{
	if (smoothtype == 0)
	{
		vtkSmoothPolyDataFilter *smoother = vtkSmoothPolyDataFilter::New();
		 smoother->SetInputData(source);
		 smoother->SetNumberOfIterations(NumIt);
		 smoother->SetRelaxationFactor(RelaxFactor);
		 smoother->SetBoundarySmoothing(BoundarySmooth);
		 smoother->SetFeatureEdgeSmoothing(FeatureEdgeSmooth);
		 smoother->SetFeatureAngle(FeatureAngle);
		 smoother->SetGenerateErrorScalars(GenerateErrScal);

		 smoother->Update();

		smooth->DeepCopy(smoother->GetOutput());
		smoother->Delete();
	}
	else if (smoothtype == 1)
	{
		vtkWindowedSincPolyDataFilter *smoother = vtkWindowedSincPolyDataFilter::New();
		smoother->SetInputData(source);
		smoother->SetNumberOfIterations(NumIt);
		smoother->SetPassBand(RelaxFactor);
		smoother->SetBoundarySmoothing(BoundarySmooth);
		smoother->SetFeatureEdgeSmoothing(FeatureEdgeSmooth);
		smoother->SetFeatureAngle(FeatureAngle);
		smoother->SetGenerateErrorScalars(GenerateErrScal);
		smoother->Update();

		smooth->DeepCopy(smoother->GetOutput());
		smoother->Delete();
	}
	else
	{
		vtkSmoothPolyDataFilter* smoother = vtkSmoothPolyDataFilter::New();
		smoother->SetInputData(source);
		smoother->SetNumberOfIterations(NumIt);
		smoother->SetRelaxationFactor(RelaxFactor);
		smoother->SetBoundarySmoothing(BoundarySmooth);
		smoother->SetFeatureEdgeSmoothing(FeatureEdgeSmooth);
		smoother->SetFeatureAngle(FeatureAngle);
		smoother->SetGenerateErrorScalars(GenerateErrScal);

		smoother->Update();

		smooth->DeepCopy(smoother->GetOutput());
		smoother->Delete();
	}
}

void CProcess3DData::DoComputeFeatureEdges(vtkPolyData *source, vtkPolyData *FeatureEdges, 
	bool boundary, bool nonManifold, bool manifold, bool sharp, double FeatureAngle)
{

	vtkFeatureEdges *feature = vtkFeatureEdges::New();
	feature->SetInputData(source);
	feature->SetBoundaryEdges(boundary);
	feature->SetManifoldEdges(manifold);
	feature->SetNonManifoldEdges(nonManifold);
	feature->SetFeatureEdges(sharp);
	feature->SetFeatureAngle(FeatureAngle);
	feature->ColoringOn();

	//if (EdgeType == 0)
	//	feature->SetBoundaryEdges(true);
	//else if (EdgeType == 1)
	//	feature->SetNonManifoldEdges(true);
	//else if (EdgeType == 2)
	//	feature->SetManifoldEdges(true);
	//else if (EdgeType == 3)
	//{
	//	feature->SetFeatureEdges(true);
	//	
	//}
	
	feature->Update();
	FeatureEdges->DeepCopy(feature->GetOutput());
	feature->Delete();
}

void CProcess3DData::ExtractOuterSurface(vtkPolyData *source, vtkPolyData *outsurf)
{
	vtkPolyDataConnectivityFilter *connect = vtkPolyDataConnectivityFilter::New();
	connect->SetInputData(source);
	connect->ColorRegionsOff();
	connect->SetExtractionModeToSpecifiedRegions();
	connect->Update();

	int N = connect->GetNumberOfExtractedRegions();

	int outsurfID = 0;
	double maxMean = 0;
		
	// Find outer surface
	for (int i = 0; i < N; i++)
	{
		connect->InitializeSpecifiedRegionList();
		connect->AddSpecifiedRegion(i);
		connect->Update();

		vtkRemoveUnusedPolyDataPoints *rem = vtkRemoveUnusedPolyDataPoints::New();
		rem->SetInputConnection(connect->GetOutputPort());
		rem->Update();

		double CMmean = 0;
		double CMSdev = 0;


		vtkExtMisc::DistanceFromCMStats(rem->GetOutput()->GetPoints(), CMmean, CMSdev);
		const int N1 = rem->GetOutput()->GetNumberOfPoints();
		const int N2 = rem->GetOutput()->GetNumberOfPolys();


		if (CMmean > maxMean)
		{
			maxMean = CMmean;
			outsurfID = i;
		}
		rem->Delete();
	}
	connect->InitializeSpecifiedRegionList();
	connect->AddSpecifiedRegion(outsurfID);
	connect->Update();

	vtkRemoveUnusedPolyDataPoints *rem = vtkRemoveUnusedPolyDataPoints::New();
	rem->SetInputConnection(connect->GetOutputPort());
	rem->Update();

	outsurf->DeepCopy(rem->GetOutput());
	connect->Delete();
	rem->Delete();
}
