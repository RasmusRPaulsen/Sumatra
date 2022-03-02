/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkConstrainedICPTransform.cxx,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkConstrainedICPTransform.h"

#include "vtkCellLocator.h"
#include "vtkPointLocator.h"
#include "vtkDataSet.h"
#include "vtkLandmarkTransform.h"
#include "vtkLandmarkTransformDirectQuart.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkTransform.h"
#include "vtkPolyData.h"
#include "vtkFeatureEdges.h"
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
//#include "multitimer.h"
#include "MeshMeasures.h"

vtkStandardNewMacro(vtkConstrainedICPTransform);

//----------------------------------------------------------------------------

vtkConstrainedICPTransform::vtkConstrainedICPTransform()
: vtkLinearTransform()
{
	featureLocator = NULL;
	featureEdges = NULL;
	this->Source = NULL;
	this->Target = NULL;
//	this->Locator = NULL;
	this->LandmarkTransform = vtkLandmarkTransform::New();
	this->MaximumNumberOfIterations = 50;
	this->CheckMeanDistance = 0;
	this->MeanDistanceMode = VTK_ICP_MODE_RMS;
	this->MaximumMeanDistance = 0.01;
	this->MaximumNumberOfLandmarks = 200;
	this->StartByMatchingCentroids = 0;

	this->NumberOfIterations = 0;
	this->MeanDistance = 0.0;
	this->LocatorType = 1;
	DebugOutput = false;

	Quarternion[0] = 0;
	Quarternion[1] = 0;
	Quarternion[2] = 0;
	Quarternion[3] = 0;
	Translation[0] = 0;
	Translation[1] = 0;
	Translation[2] = 0;
	CheckMatrix = vtkMatrix4x4::New();

	EdgeFinderTolerance = 0.000000001;
	CheckDistance = 0;
	DistanceThreshold = 10;
	MeanDistanceAfterThreshold = 0;
}

//----------------------------------------------------------------------------

const char *vtkConstrainedICPTransform::GetMeanDistanceModeAsString()
{
	if ( this->MeanDistanceMode == VTK_ICP_MODE_RMS )
	{
		return "RMS";
	}
	else
	{
		return "AbsoluteValue";
	}
}

//----------------------------------------------------------------------------

vtkConstrainedICPTransform::~vtkConstrainedICPTransform()
{
	ReleaseSource();
	ReleaseTarget();
// 	ReleaseLocator();
	this->LandmarkTransform->Delete();
	if (featureEdges)
		featureEdges->Delete();
	if (featureLocator)
		featureLocator->Delete();
	CheckMatrix->Delete();
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::SetSource(vtkDataSet *source)
{
	if (this->Source == source)
	{
		return;
	}

	if (this->Source)
	{
		this->ReleaseSource();
	}

	if (source)
	{
		source->Register(this);
	}

	this->Source = source;
	this->Modified();
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::ReleaseSource(void) {
	if (this->Source) 
	{
		this->Source->UnRegister(this);
		this->Source = NULL;
	}
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::SetTarget(vtkDataSet *target)
{
	if (this->Target == target)
	{
		return;
	}

	if (this->Target)
	{
		this->ReleaseTarget();
	}

	if (target)
	{
		target->Register(this);
	}

	this->Target = target;
	this->Modified();
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::ReleaseTarget(void) {
	if (this->Target) 
	{
		this->Target->UnRegister(this);
		this->Target = NULL;
	}
}

//----------------------------------------------------------------------------

//void vtkConstrainedICPTransform::SetLocator(vtkCellLocator *locator)
//{
//	if (this->Locator == locator)
//	{
//		return;
//	}
//
//	if (this->Locator)
//	{
//		this->ReleaseLocator();
//	}
//
//	if (locator)
//	{
//		locator->Register(this);
//	}
//
//	this->Locator = locator;
//	this->Modified();
//}

//----------------------------------------------------------------------------

//void vtkConstrainedICPTransform::ReleaseLocator(void) {
//	if (this->Locator) 
//	{
//	this->Locator->UnRegister(this);
//		this->Locator = NULL;
//	}
//}

//----------------------------------------------------------------------------

//void vtkConstrainedICPTransform::CreateDefaultLocator() {
//	if (this->Locator) 
//	{
//		this->ReleaseLocator();
//	}
//
//	this->Locator = vtkCellLocator::New();
//}

//------------------------------------------------------------------------

vtkMTimeType vtkConstrainedICPTransform::GetMTime()
{
	unsigned long result = this->vtkLinearTransform::GetMTime();
	unsigned long mtime;

	if (this->Source)
	{
		mtime = this->Source->GetMTime(); 
		if (mtime > result)
		{
			result = mtime;
		}
	}

	if (this->Target)
	{
		mtime = this->Target->GetMTime(); 
		if (mtime > result)
		{
			result = mtime;
		}
	}

	//if (this->Locator)
	//{
	//	mtime = this->Locator->GetMTime(); 
	//	if (mtime > result)
	//	{
	//		result = mtime;
	//	}
	//}

	if (this->LandmarkTransform)
	{
		mtime = this->LandmarkTransform->GetMTime();
		if (mtime > result)
		{
			result = mtime;
		}
	}

	return result;
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::Inverse()
{
	vtkDataSet *tmp1 = this->Source;
	this->Source = this->Target;
	this->Target = tmp1;
	this->Modified();
}

//----------------------------------------------------------------------------

vtkAbstractTransform *vtkConstrainedICPTransform::MakeTransform()
{
	return vtkConstrainedICPTransform::New(); 
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::InternalDeepCopy(vtkAbstractTransform *transform)
{
	vtkConstrainedICPTransform *t = (vtkConstrainedICPTransform *)transform;

	this->SetSource(t->GetSource());
	this->SetTarget(t->GetTarget());
//	this->SetLocator(t->GetLocator());
	this->SetMaximumNumberOfIterations(t->GetMaximumNumberOfIterations());
	this->SetCheckMeanDistance(t->GetCheckMeanDistance());
	this->SetMeanDistanceMode(t->GetMeanDistanceMode());
	this->SetMaximumMeanDistance(t->GetMaximumMeanDistance());
	this->SetMaximumNumberOfLandmarks(t->GetMaximumNumberOfLandmarks());

	this->Modified();
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::InternalUpdate()
{
	// Check source, target

	if (this->Source == NULL || !this->Source->GetNumberOfPoints())
	{
		vtkErrorMacro(<<"Can't execute with NULL or empty input");
		return;
	}

	if (this->Target == NULL || !this->Target->GetNumberOfPoints())
	{
		vtkErrorMacro(<<"Can't execute with NULL or empty target");
		return;
	}

//	CMultiTimer::MultiTimer().Start("ICP CreateFeatureEdges", 1);
	CreateFeatureEdges();
//	CMultiTimer::MultiTimer().End("ICP CreateFeatureEdges");

	// Create locator
//	this->CreateDefaultLocator();

	vtkCellLocator *CellLocator = NULL;
	vtkPointLocator *PointLocator = NULL;

//	CMultiTimer::MultiTimer().Start("Build Locator", 1);
	
	// Do they last few iterations with point to surface mode (if in point to point mode)
	bool FinishWithPointToSurface = true;
	if (LocatorType == 0)
		FinishWithPointToSurface = false;

	if (LocatorType == 0 || FinishWithPointToSurface)
	{
		CellLocator = vtkCellLocator::New();
		CellLocator->SetDataSet(this->Target);
		CellLocator->SetNumberOfCellsPerBucket(1);
		CellLocator->BuildLocator();
	}

	if (LocatorType == 1)
	{
		PointLocator = vtkPointLocator::New();
		PointLocator->SetDataSet(this->Target);
		PointLocator->SetNumberOfPointsPerBucket(1);
		PointLocator->BuildLocator();
	}

//	CMultiTimer::MultiTimer().End("Build Locator");
	// Create two sets of points to handle iteration

	int step = 1;
	if (this->Source->GetNumberOfPoints() > this->MaximumNumberOfLandmarks)
	{
		step = this->Source->GetNumberOfPoints() / this->MaximumNumberOfLandmarks;
		vtkDebugMacro(<< "Landmarks step is now : " << step);
	}

	vtkIdType nb_points = this->Source->GetNumberOfPoints() / step;

	// Allocate some points. ... Rasmus hack: Not really valid any more
	// - closestp is used so that the internal state of LandmarkTransform remains
	//   correct whenever the iteration process is stopped (hence its source
	//   and landmark points might be used in a vtkThinPlateSplineTransform).
	// - points2 could have been avoided, but do not ask me why 
	//   InternalTransformPoint is not working correctly on my computer when
	//   in and out are the same pointer.

	vtkPoints *points1 = vtkPoints::New();
	points1->SetNumberOfPoints(nb_points);

	// we want to find the final transformation as a quarternion that is not a result of a chain of concatenations
	// For hacky reasons we keep an original version of the used source points
	vtkPoints *QuartHackPointsOrg = vtkPoints::New();
	QuartHackPointsOrg->SetNumberOfPoints(nb_points);

	//vtkPoints *closestp = vtkPoints::New();
	//closestp->SetNumberOfPoints(nb_points);

	vtkPoints *points2 = vtkPoints::New();
	points2->SetNumberOfPoints(nb_points);

	// Fill with initial positions (sample dataset using step)
	vtkTransform *accumulate = vtkTransform::New();
	accumulate->PostMultiply();

	vtkIdType i;
	int j;
	double p1[3], p2[3];

	if (DebugOutput)
		std::cout << "Constrained ICP with max iterations " << this->MaximumNumberOfIterations << " Max landmarks " << this->MaximumNumberOfLandmarks << " max meanerror " << this->MaximumMeanDistance << std::endl;
	if (StartByMatchingCentroids)
	{
//		CMultiTimer::MultiTimer().Start("ICP Match Centroids", 1);
		if (DebugOutput)
			std::cout << "Starting ICP by matching centroids" << std::endl;
		double source_centroid[3] = {0,0,0};
		for (i = 0; i < this->Source->GetNumberOfPoints(); i++)
		{
			this->Source->GetPoint(i, p1);
			source_centroid[0] += p1[0];
			source_centroid[1] += p1[1];
			source_centroid[2] += p1[2];
		}
		source_centroid[0] /= this->Source->GetNumberOfPoints();
		source_centroid[1] /= this->Source->GetNumberOfPoints();
		source_centroid[2] /= this->Source->GetNumberOfPoints();

		double target_centroid[3] = {0,0,0};
		for (i = 0; i < this->Target->GetNumberOfPoints(); i++)
		{
			this->Target->GetPoint(i, p2);
			target_centroid[0] += p2[0];
			target_centroid[1] += p2[1];
			target_centroid[2] += p2[2];
		}
		target_centroid[0] /= this->Target->GetNumberOfPoints();
		target_centroid[1] /= this->Target->GetNumberOfPoints();
		target_centroid[2] /= this->Target->GetNumberOfPoints();

		accumulate->Translate(target_centroid[0] - source_centroid[0],
			target_centroid[1] - source_centroid[1],
			target_centroid[2] - source_centroid[2]);
		accumulate->Update();

		for (i = 0, j = 0; i < nb_points; i++, j += step)
		{
			double outPoint[3];
			accumulate->InternalTransformPoint(this->Source->GetPoint(j),
				outPoint);
			points1->SetPoint(i, outPoint);
		}
//		CMultiTimer::MultiTimer().End("ICP Match Centroids");
	}
	else 
	{
		if (DebugOutput)
			std::cout << "ICP not matching centroids" << std::endl;
		for (i = 0, j = 0; i < nb_points; i++, j += step)
		{
			points1->SetPoint(i, this->Source->GetPoint(j));
		}
	}

	// Get an original copy - no matter if CM is subtracted or not
	for (i = 0, j = 0; i < nb_points; i++, j += step)
	{
		QuartHackPointsOrg->SetPoint(i, this->Source->GetPoint(j));
	}


	// Go

	vtkIdType cell_id;
	int sub_id;
	double dist2, totaldist = 0;
	double outPoint[3];

	vtkPoints *temp, *a = points1, *b = points2;

	this->NumberOfIterations = 0;


	vtkPoints *QuartHackPoints1 = vtkPoints::New(); // needs to come from the original source point set
	vtkPoints *QuartHackPoints2 = vtkPoints::New(); // needs to come from the latest target point set

	int NMarksUsed = nb_points * 0.03;
	do 
	{
		if (NumberOfIterations > 0.1 * MaximumNumberOfIterations)
			NMarksUsed = nb_points * 0.05;
		if (NumberOfIterations > 0.2 * MaximumNumberOfIterations)
			NMarksUsed = nb_points * 0.1;
		if (NumberOfIterations > 0.3 * MaximumNumberOfIterations)
			NMarksUsed = nb_points;

		step = nb_points / NMarksUsed;
		vtkIdType nb_points2 = nb_points / step;

		// Fill points with the closest points to each vertex in input
		vtkPoints *pp1 = vtkPoints::New();
		vtkPoints *pp2 = vtkPoints::New();
		QuartHackPoints1->Reset();
		QuartHackPoints2->Reset();

//		CMultiTimer::MultiTimer().Start("ICP Closest Points", 1);
		int borderpoints = 0;
		int distancePoints = 0;
		int validPoints = 0;
		MeanDistanceAfterThreshold = 0;

		for(i = 0,j=0; i < nb_points2; i++, j += step)
		{
			if (LocatorType == 0)
			{
				CellLocator->FindClosestPoint(a->GetPoint(j),
					outPoint,
					cell_id,
					sub_id,
					dist2);

				bool valid = true;
				double tcp[3];
				// see if it is on a border
				featureLocator->FindClosestPoint(outPoint, tcp, cell_id, sub_id, dist2);
				if (dist2 < EdgeFinderTolerance)
				{
					borderpoints++;
					valid = false;
				}
				if (CheckDistance)
				{
					if (dist2 > DistanceThreshold*DistanceThreshold)
					{
						valid = false;
						distancePoints++;
					}
				}
				if (valid)
				{
					pp1->InsertNextPoint(a->GetPoint(j));
					pp2->InsertNextPoint(outPoint);

					QuartHackPoints1->InsertNextPoint(QuartHackPointsOrg->GetPoint(j));
					QuartHackPoints2->InsertNextPoint(outPoint);

					MeanDistanceAfterThreshold += sqrt(dist2);
					validPoints++;
				}
			}
			else
			{
				vtkIdType pid = PointLocator->FindClosestPoint(a->GetPoint(j));

				Target->GetPoint(pid, outPoint);
				// see if it is on a border
				if (BorderPointsID[pid])
				{
					borderpoints++;
				}
				else
				{
					pp1->InsertNextPoint(a->GetPoint(j));
					pp2->InsertNextPoint(outPoint);

					QuartHackPoints1->InsertNextPoint(QuartHackPointsOrg->GetPoint(j));
					QuartHackPoints2->InsertNextPoint(outPoint);
				}

			}
		}
//		CMultiTimer::MultiTimer().End("ICP Closest Points");

		if (DebugOutput)
			std::cout << "# " << NumberOfIterations << " Points: " << nb_points2 << " Borderpoints: " << borderpoints << " non-borders: " << pp1->GetNumberOfPoints() << " distpoints " << distancePoints << std::endl;

		vtkDebugMacro(<< "Borderpoints: " << borderpoints << " non-borders: " << pp1->GetNumberOfPoints());

		if (validPoints)
		{
			MeanDistanceAfterThreshold /= validPoints;
		}
		if (DebugOutput)
		{
			std::cout << "Meand distance after threshold " << MeanDistanceAfterThreshold << std::endl;
		}

//		CMultiTimer::MultiTimer().Start("ICP Update Transform", 1);
		// Build the landmark transform
		this->LandmarkTransform->SetSourceLandmarks(pp1);
		this->LandmarkTransform->SetTargetLandmarks(pp2);
		this->LandmarkTransform->Update();

		// Concatenate (can't use this->Concatenate directly)
		accumulate->Concatenate(this->LandmarkTransform->GetMatrix());
//		CMultiTimer::MultiTimer().End("ICP Update Transform");

		this->NumberOfIterations++;
		vtkDebugMacro(<< "Iteration: " << this->NumberOfIterations);
		if (this->NumberOfIterations >= this->MaximumNumberOfIterations) 
		{
			if (FinishWithPointToSurface && LocatorType == 1)
			{
				if (DebugOutput)
					std::cout << "Switching to point to surface mode" << std::endl;
				LocatorType = 0;
			}
			else
			{
				if (DebugOutput)
					std::cerr << "ICP stopped since max iterations reached: " << this->MaximumNumberOfIterations << std::endl;
				break;
			}
		}

		// Move mesh and compute mean distance if needed
		if (this->CheckMeanDistance)
		{
			totaldist = 0.0;
		}

		
//		CMultiTimer::MultiTimer().Start("ICP Move Mesh", 1);
		for(i = 0; i < nb_points; i++)
		{
			a->GetPoint(i, p1);
			this->LandmarkTransform->InternalTransformPoint(p1, p2);
			b->SetPoint(i, p2);
			if (this->CheckMeanDistance)
			{
				if (this->MeanDistanceMode == VTK_ICP_MODE_RMS) 
				{
					totaldist += vtkMath::Distance2BetweenPoints(p1, p2);
				} else {
					totaldist += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
				}
			}
		}
//		CMultiTimer::MultiTimer().End("ICP Move Mesh");

		// We can only step if the full number of landmarks has been tried at least once
		if (this->CheckMeanDistance && (NMarksUsed == nb_points))
		{
			if (this->MeanDistanceMode == VTK_ICP_MODE_RMS) 
			{
				this->MeanDistance = sqrt(totaldist / (double)nb_points);
			} else {
				this->MeanDistance = totaldist / (double)nb_points;
			}
			vtkDebugMacro("Mean distance: " << this->MeanDistance);
			if (DebugOutput)
				std::cout << "Mean distance: " << this->MeanDistance << std::endl;
	
			if (this->MeanDistance <= this->MaximumMeanDistance)
			{
				if (FinishWithPointToSurface && LocatorType == 1)
				{
					if (DebugOutput)
						std::cout << "Switching to point to surface mode" << std::endl;
					LocatorType = 0;
				}
				else
				{
					if (DebugOutput)
						std::cerr << "ICP stopped since mean distance reached: " << this->MaximumMeanDistance << std::endl;
					break;
				}
			}
		}

		temp = a;
		a = b;
		b = temp;
	
		pp1->Delete();
		pp2->Delete();
	} 
	while (1);

	// Now recover accumulated result

	this->Matrix->DeepCopy(accumulate->GetMatrix());

//	CMultiTimer::MultiTimer().Start("ICP Compute Final Transform", 1);
	vtkLandmarkTransformDirectQuart *QuartTrans = vtkLandmarkTransformDirectQuart::New();
	QuartTrans->SetSourceLandmarks(QuartHackPoints1);
	QuartTrans->SetTargetLandmarks(QuartHackPoints2);
	QuartTrans->Update();

	CheckMatrix->DeepCopy(QuartTrans->GetMatrix());
		
	for (int i = 0; i < 3; i++)
	{
		Translation[i] = QuartTrans->GetTranslation()[i];
	}
	for (int i = 0; i < 4; i++)
	{
		Quarternion[i] = QuartTrans->GetQuarternion()[i];
	}
	QuartTrans->Delete();
//	CMultiTimer::MultiTimer().End("ICP Compute Final Transform");

//	std::cout << "Quartenion hack: set 1 " << QuartHackPoints1->GetNumberOfPoints() << " set 2 " << QuartHackPoints2->GetNumberOfPoints() << std::endl;


	if (CellLocator)
	{
		CellLocator->Delete();
	}
	if (PointLocator)
	{
		PointLocator->Delete();
	}
	accumulate->Delete();
	points1->Delete();
//	closestp->Delete();
	points2->Delete();

	QuartHackPointsOrg->Delete();
	QuartHackPoints1->Delete();
	QuartHackPoints2->Delete();
}

//----------------------------------------------------------------------------

void vtkConstrainedICPTransform::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);

	if ( this->Source ) 
	{
		os << indent << "Source: " << this->Source << "\n";
	}
	else 
	{
		os << indent << "Source: (none)\n";
	}

	if ( this->Target ) 
	{
		os << indent << "Target: " << this->Target << "\n";
	}
	else 
	{
		os << indent << "Target: (none)\n";
	}

	//if ( this->Locator ) 
	//{
	//	os << indent << "Locator: " << this->Locator << "\n";
	//}
	//else 
	//{
	//	os << indent << "Locator: (none)\n";
	//}

	os << indent << "MaximumNumberOfIterations: " << this->MaximumNumberOfIterations << "\n";
	os << indent << "CheckMeanDistance: " << this->CheckMeanDistance << "\n";
	os << indent << "MeanDistanceMode: " << this->GetMeanDistanceModeAsString() << "\n";
	os << indent << "MaximumMeanDistance: " << this->MaximumMeanDistance << "\n";
	os << indent << "MaximumNumberOfLandmarks: " << this->MaximumNumberOfLandmarks << "\n";
	os << indent << "StartByMatchingCentroids: " << this->StartByMatchingCentroids << "\n";
	os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";
	os << indent << "MeanDistance: " << this->MeanDistance << "\n";
	if(this->LandmarkTransform)
	{
		os << indent << "LandmarkTransform:\n";
		this->LandmarkTransform->PrintSelf(os, indent.GetNextIndent());
	}
}

void vtkConstrainedICPTransform::CreateFeatureEdges()
{
	vtkDebugMacro(<< "CreateFeatureEdges");

	vtkFeatureEdges *feature = vtkFeatureEdges::New();
	feature->SetInputData(Target);
	feature->BoundaryEdgesOn();
	feature->NonManifoldEdgesOff();
	feature->FeatureEdgesOff();
	feature->Update();

	if (DebugOutput)
		std::cout << "Border edges " << feature->GetOutput()->GetNumberOfLines() << std::endl;

	featureEdges = vtkPolyData::New();
	featureEdges->DeepCopy(feature->GetOutput());

	featureLocator = vtkCellLocator::New();
	featureLocator->SetDataSet(featureEdges);
	featureLocator->SetNumberOfCellsPerBucket(1);
	featureLocator->BuildLocator();

	feature->Delete();

	// Compute the EdgeFinder tolerance to be half the median interpoint distance
	double EdgeFinderTolerance = 0.00000001;

	if (DebugOutput)
		std::cout << "Mesh statistics" << std::endl;
	double minD = 0;
	double maxD = 0;
	double meanD = 0;
	double sdvD = 0;
	double medD = 0;
	double D05 = 0;
	double D95 = 0;

	CMeshMeasures::NearestNeighbourStatistics(Target, minD, maxD, meanD, sdvD, medD, D05, D95);

	EdgeFinderTolerance = medD / 2;
	
	if (DebugOutput)
		std::cout << "Edgefinder distance tolerance " << EdgeFinderTolerance << std::endl;

	// Mark points if they are on the border
	// very hacky way!
	BorderPointsID.resize(Target->GetNumberOfPoints());

	vtkIdType cell_id;
	int sub_id;
	double dist2, totaldist = 0;
	double P[3];

	for (int i = 0; i < Target->GetNumberOfPoints(); i++)
	{
		Target->GetPoint(i, P);

		double tcp[3];
		// see if it is on a border
		featureLocator->FindClosestPoint(P, tcp, cell_id, sub_id, dist2);
		if (dist2 < EdgeFinderTolerance)
		{
			BorderPointsID[i] = 1;
		}
		else
		{
			BorderPointsID[i] = 0;
		}
	}
}
