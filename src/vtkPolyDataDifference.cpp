#include "vtkPolyDataDifference.h"


#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkDoubleArray.h"
#include <vtkTriangle.h>
#include <vtkFeatureEdges.h>
#include <iostream>
#include "MeshMeasures.h"

#include <algorithm>


vtkStandardNewMacro(vtkPolyDataDifference);

//----------------------------------------------------------------------------
vtkPolyDataDifference::vtkPolyDataDifference()
{
	featureLocator = NULL;
	featureEdges = NULL;

	SignedDistance = false;
	ExcludeEdgePoints = false;
	MaximumNumberOfLandmarks = -1;
	EdgeFinderTolerance = 0.000000001;
	this->SetNumberOfInputPorts(2);

}
//----------------------------------------------------------------------------
vtkPolyDataDifference::~vtkPolyDataDifference()
{
	if (featureEdges)
		featureEdges->Delete();
	if (featureLocator)
		featureLocator->Delete();
}

void vtkPolyDataDifference::SetTargetData(vtkPolyData *target)
{
  this->SetInputData(1, target);
}

vtkPolyData *vtkPolyDataDifference::GetTarget()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkPolyData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
}


////----------------------------------------------------------------------------
//void vtkPolyDataDifference::SetTarget(vtkPolyData *pd)
//{
//  this->vtkProcessObject::SetNthInput(1,pd);
//}
//
////----------------------------------------------------------------------------
//vtkPolyData *vtkPolyDataDifference::GetTarget()
//{
//  if (this->NumberOfInputs < 2) {
//	return NULL;
//  }
//  return (vtkPolyData *)(this->Inputs[1]);
//}


void vtkPolyDataDifference::DifferenceWithEdgePointsExcluded(vtkInformationVector **inputVector,  vtkInformationVector *outputVector)
{
//	std::cout << "Difference with edge points excluded" << std::endl;

    // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *targetInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  vtkPolyData *target = 0;
  if (targetInfo)
    {
    target = vtkPolyData::SafeDownCast(
      targetInfo->Get(vtkDataObject::DATA_OBJECT()));
    }


	CreateFeatureEdges(target);

	vtkCellLocator *locator = vtkCellLocator::New();

	locator->SetDataSet(target);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	//vtkDoubleArray* newnorms = vtkDoubleArray::New();
	//newnorms->SetNumberOfComponents(3);
	//newnorms->SetName("Normals");

	vtkDataArray *normals = target->GetPointData()->GetNormals();

	int NPoints = input->GetPoints()->GetNumberOfPoints();

	if (MaximumNumberOfLandmarks > -1)
	{
		NPoints = std::min(MaximumNumberOfLandmarks, NPoints);
	}
	int step = input->GetNumberOfPoints() / NPoints;

	int borderpoints = 0;
	int count = 0;
	for (int k2 = 0, k = 0; k2 < NPoints; k2++, k += step)
	{
//		std::cout << input->GetNumberOfPoints() << " k " << k << " k2 " << k2 << std::endl;
		int id = k; 

		double P[3];
		input->GetPoint(id, P);

		// closest point
		double cp[3];

		int sub_id;
		double dist2 = 0;
		vtkIdType cell_id;

		locator->FindClosestPoint(P,
			cp,
			cell_id,
			sub_id,
			dist2);

		double td = sqrt(dist2);

		if (SignedDistance && normals != NULL)
		{
			// Find weighted normal
			vtkTriangle *tri = vtkTriangle::SafeDownCast(target->GetCell(cell_id));
			if (!tri)
			{
				vtkErrorMacro("Update: Only able to handle triangles");
			}
			else
			{
				double norm[3];

				// Closest point in evaluation
				double p2[3];

				// Parametric coordinates
				double pcoords[3];

				// Distance
				double dist2;

				double weights[3];

				int inside = tri->EvaluatePosition(P,p2,sub_id,pcoords,dist2,weights);

				int id0 = tri->GetPointId(0);
				int id1 = tri->GetPointId(1);
				int id2 = tri->GetPointId(2);

				double norm0[3];
				normals->GetTuple(id0, norm0);
				double norm1[3];
				normals->GetTuple(id1, norm1);
				double norm2[3];
				normals->GetTuple(id2, norm2);

				// Find x coordinate by weighting x coordinates off the three normals (doing the same with y and z)
				norm[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
				norm[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
				norm[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];

				double vec[3];
				vec[0] = P[0] - cp[0];
				vec[1] = P[1] - cp[1];
				vec[2] = P[2] - cp[2];

				double dot1 = (norm[0] * vec[0] + norm[1] * vec[1] + norm[2] * vec[2]);

				if (dot1 < 0)
					td = -td;
			}
		}

		double tcp[3];
		// see if it is on a border
		featureLocator->FindClosestPoint(cp, tcp, cell_id, sub_id, dist2);
		if (dist2 < EdgeFinderTolerance)
		{
			borderpoints++;
		}
		else
		{
			pts->InsertNextPoint(P);
			int tid2 = verts->InsertNextCell(1);
			verts->InsertCellPoint(tid2);
			scalars->InsertNextTuple(&td);

			//if (normals)
			//{
			//	double n[3];
			//	normals->GetTuple(k, n);
			//	newnorms->InsertNextTuple(n);
			//}

		}
	}
//	std::cout << "Difference: Borderpoints: " << borderpoints << " non-borders: " << pts->GetNumberOfPoints() << std::endl;

	output->SetPoints(pts);
	output->SetVerts(verts);
	output->GetPointData()->SetScalars(scalars);
	output->GetPointData()->SetNormals(NULL);
	//if (normals)
	//	output->GetPointData()->SetNormals(newnorms);

	scalars->Delete();
	//newnorms->Delete();
	pts->Delete();
	verts->Delete();
	locator->Delete();
}

//----------------------------------------------------------------------------
int vtkPolyDataDifference::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  if (ExcludeEdgePoints)
	{
		DifferenceWithEdgePointsExcluded(inputVector, outputVector);
		return 1;
	}

    // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *targetInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  vtkPolyData *target = 0;
  if (targetInfo)
    {
    target = vtkPolyData::SafeDownCast(
      targetInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

//	vtkDebugMacro(<< "-->Starting vtkPolyDataDifference");
	
	//vtkPolyData *input		 = this->GetInput();
	//vtkPolyData *output		 = this->GetOutput();
	//vtkPolyData *target		 = this->GetTarget();
	
	output->CopyStructure(input);
	output->GetPointData()->PassData(input->GetPointData());
	output->GetCellData()->PassData(input->GetCellData());

	vtkCellLocator *locator = vtkCellLocator::New();

	locator->SetDataSet(target);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();


	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

 	vtkDataArray *normals = target->GetPointData()->GetNormals();

	int count = 0;
	for (int k = 0; k < input->GetPoints()->GetNumberOfPoints(); k++)
	{
//		if (!(count++ % 500))
//			std::cout << k << " / " << input->GetPoints()->GetNumberOfPoints() << std::endl;

		int id = k; 

		double P[3];
		input->GetPoint(id, P);

		// closest point
		double cp[3];
			
		int sub_id;
		double dist2 = 0;
		vtkIdType cell_id;

		locator->FindClosestPoint(P,
									cp,
									cell_id,
									sub_id,
									dist2);
		
		

		double td = sqrt(dist2);

		if (SignedDistance && normals != NULL)
		{
			// Find weighted normal
			vtkTriangle *tri = vtkTriangle::SafeDownCast(target->GetCell(cell_id));
			if (!tri)
			{
			   vtkErrorMacro("Update: Only able to handle triangles");
			}
			else
			{
				double norm[3];

				// Closest point in evaluation
				double p2[3];
				
				// Parametric coordinates
				double pcoords[3];
				
				// Distance
				double dist2;
				
				double weights[3];
				
				int inside = tri->EvaluatePosition(P,p2,sub_id,pcoords,dist2,weights);
				
				int id0 = tri->GetPointId(0);
				int id1 = tri->GetPointId(1);
				int id2 = tri->GetPointId(2);
				
				double norm0[3];
				normals->GetTuple(id0, norm0);
				double norm1[3];
				normals->GetTuple(id1, norm1);
				double norm2[3];
				normals->GetTuple(id2, norm2);
				
				// Find x coordinate by weighting x coordinates off the three normals (doing the same with y and z)
				norm[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
				norm[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
				norm[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];

				double vec[3];
				vec[0] = P[0] - cp[0];
				vec[1] = P[1] - cp[1];
				vec[2] = P[2] - cp[2];

				double dot1 = (norm[0] * vec[0] + norm[1] * vec[1] + norm[2] * vec[2]);

				if (dot1 < 0)
					td = -td;
			}
		}

		scalars->InsertTuple(id, &td);
	}


	output->GetPointData()->SetScalars(scalars);
	scalars->Delete();
	vtkDebugMacro(<< "<--Ending vtkPolyDataDifference");
	return 1;
}

/*
void CVectorFieldDiffusion::FindWeigthedNormal(double *tp,double *normal)
{
	vtkIdType cellId;
	// Closest point on surface
	double cp[3];
	double dist;
	int subId;
	Locator->FindClosestPoint(tp, cp, cellId, subId, dist);
	
	vtkTriangle *tri = vtkTriangle::SafeDownCast(target->GetCell(cellId));
	if (!tri)
	{
		//    vtkErrorMacro("Update: Only able to handle triangles");
		return;
	}
	
	// Closest point in evaluation
	double p2[3];
	
	// Parametric coordinates
	double pcoords[3];
	
	// Distance
	double dist2;
	
	double weights[3];
	
	int inside = tri->EvaluatePosition(tp,p2,subId,pcoords,dist2,weights);
//	std::cout << "Weights: " << weights[0] << ", " << weights[1] << ", " << weights[2] << std::endl;
	
	int id0 = tri->GetPointId(0);
	int id1 = tri->GetPointId(1);
	int id2 = tri->GetPointId(2);
	
//	std::cout << "id: " << id0 << ", " << id1 << ", " << id2 << std::endl;
	
	vtkDataArray *normals = target->GetPointData()->GetNormals(); 
	double norm0[3];
	normals->GetTuple(id0, norm0);
	double norm1[3];
	normals->GetTuple(id1, norm1);
	double norm2[3];
	normals->GetTuple(id2, norm2);
	
	// Find x coordinate by weihgting x coordinates off the three normals (doing the same with y and z)
	normal[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
	normal[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
	normal[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];
	
}
*/

void vtkPolyDataDifference::CreateFeatureEdges(vtkPolyData *target)
{
	vtkDebugMacro(<< "CreateFeatureEdges");

	vtkFeatureEdges *feature = vtkFeatureEdges::New();
	feature->SetInputData(target);
	feature->BoundaryEdgesOn();
	feature->NonManifoldEdgesOff();
	feature->FeatureEdgesOff();
	feature->Update();

//	std::cout << "Border edges " << feature->GetOutput()->GetNumberOfLines() << std::endl;

	featureEdges = vtkPolyData::New();
	featureEdges->DeepCopy(feature->GetOutput());

	featureLocator = vtkCellLocator::New();
	featureLocator->SetDataSet(featureEdges);
	featureLocator->SetNumberOfCellsPerBucket(1);
	featureLocator->BuildLocator();

	feature->Delete();

	// Compute the EdgeFinder tolerance to be half the median interpoint distance
	double EdgeFinderTolerance = 0.00000001;

	double minD = 0;
	double maxD = 0;
	double meanD = 0;
	double sdvD = 0;
	double medD = 0;
	double D05 = 0;
	double D95 = 0;

	CMeshMeasures::NearestNeighbourStatistics(target, minD, maxD, meanD, sdvD, medD, D05, D95);

	EdgeFinderTolerance = medD / 2;
}

int vtkPolyDataDifference::FillInputPortInformation(int port,
                                                      vtkInformation *info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }

  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  return 1;
}


void vtkPolyDataDifference::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
