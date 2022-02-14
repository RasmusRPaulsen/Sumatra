#include "vtkPolyDataRandomResamplePoint.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include <vtkTriangle.h>
#include <iostream>
#include <vtkMath.h>

vtkStandardNewMacro(vtkPolyDataRandomResamplePoint);

//----------------------------------------------------------------------------
vtkPolyDataRandomResamplePoint::vtkPolyDataRandomResamplePoint()
{
	ResamplePoints = 1000;
	AddBorder = 5;
	CreatePointNormals = 0;
}
//----------------------------------------------------------------------------
vtkPolyDataRandomResamplePoint::~vtkPolyDataRandomResamplePoint()
{

}



//----------------------------------------------------------------------------
int vtkPolyDataRandomResamplePoint::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
    // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

	int i;

//	output->CopyStructure(input);
//	output->GetPointData()->PassData(input->GetPointData());
//	output->GetCellData()->PassData(input->GetCellData());

	double bounds[6];
	for(i=0;i<3;i++)
	{
		bounds[i*2]=input->GetBounds()[i*2];
		bounds[i*2+1]=input->GetBounds()[i*2+1];
	}

	
	vtkCellLocator *locator = vtkCellLocator::New();

	locator->SetDataSet(input);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();

	vtkDataArray *norms = input->GetPointData()->GetNormals(); 

	vtkDoubleArray *Normals = NULL;
	if (CreatePointNormals && norms)
	{
		Normals = vtkDoubleArray::New();
		Normals->SetNumberOfComponents(3);
		Normals->SetName("Normals");
	}

	
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();

	for (int k = 0; k < ResamplePoints; k++)
	{

		double P[3];
		P[0] = vtkMath::Random(bounds[0]-AddBorder, bounds[1]+AddBorder);
		P[1] = vtkMath::Random(bounds[2]-AddBorder, bounds[3]+AddBorder);
		P[2] = vtkMath::Random(bounds[4]-AddBorder, bounds[5]+AddBorder);

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
		
		vtkIdType id = pts->InsertNextPoint(cp);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);

		if (Normals)
		{
			// Find weighted normal
			vtkTriangle *tri = vtkTriangle::SafeDownCast(input->GetCell(cell_id));
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
				norms->GetTuple(id0, norm0);
				double norm1[3];
				norms->GetTuple(id1, norm1);
				double norm2[3];
				norms->GetTuple(id2, norm2);

				// Find x coordinate by weighting x coordinates off the three normals (doing the same with y and z)
				norm[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
				norm[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
				norm[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];

				Normals->InsertNextTuple(norm);
			}
		}

	}

	output->SetPoints(pts);
	pts->Delete();
	output->SetVerts(verts);
	verts->Delete();

	if (Normals)
	{
		output->GetPointData()->SetNormals(Normals);
		Normals->Delete();
	}
	
	vtkDebugMacro(<< "<--Ending vtkPolyDataRandomResamplePoint");
	return 1;
}

void vtkPolyDataRandomResamplePoint::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
