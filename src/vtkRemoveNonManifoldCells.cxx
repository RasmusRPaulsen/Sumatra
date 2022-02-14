#include "vtkRemoveNonManifoldCells.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkMergePoints.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"

vtkStandardNewMacro(vtkRemoveNonManifoldCells);

vtkRemoveNonManifoldCells::vtkRemoveNonManifoldCells()
{
}

//--------------------------------------------------------------------------
vtkRemoveNonManifoldCells::~vtkRemoveNonManifoldCells()
{
}


//--------------------------------------------------------------------------
int vtkRemoveNonManifoldCells::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkPoints   *inPts = input->GetPoints();
	vtkIdType   numPts = input->GetNumberOfPoints();

	vtkDebugMacro(<<"Beginning vtkRemoveNonManifoldCells");
	if ( (numPts<1) || (inPts == NULL ) )
	{
		vtkDebugMacro(<<"No data to Operate On!");
		return 1;
	}

	input->BuildLinks();

	vtkPoints *pts = vtkPoints::New();
	vtkDoubleArray *Normals = NULL;

	pts->SetNumberOfPoints(numPts);

	// Brutal copy
	for (int i = 0; i < numPts; i++)
	{
		double p[3];
		input->GetPoint(i, p);
		pts->SetPoint(i, p);
	}

	int ntriangles = input->GetNumberOfCells();

	vtkCellArray *newPolys = NULL;
	newPolys = vtkCellArray::New();

	vtkIdList *neighbors = vtkIdList::New();
	neighbors->Allocate(VTK_CELL_SIZE);

	int NnonManifolds = 0;

	// Only copy triangles that are not connected to a non-manifold
	for (int f = 0; f < ntriangles; f++)
	{
		vtkTriangle *tri = vtkTriangle::SafeDownCast(input->GetCell(f));
		if (tri == NULL)
		{
			vtkErrorMacro(<<"Only handles triangles");
			return 1;
		}

		int id0 = tri->GetPointId(0);
		int id1 = tri->GetPointId(1);
		int id2 = tri->GetPointId(2);

		// First check if it is a non-manifold triangle
		input->GetCellEdgeNeighbors(f,id0,id1, neighbors);
		int numNei1 = neighbors->GetNumberOfIds();
		input->GetCellEdgeNeighbors(f,id1,id2, neighbors);
		int numNei2 = neighbors->GetNumberOfIds();
		input->GetCellEdgeNeighbors(f,id2,id0, neighbors);
		int numNei3 = neighbors->GetNumberOfIds();

		if (numNei1 < 2 && numNei2 < 2 && numNei3 < 2)
		{
			newPolys->InsertNextCell(3);

			newPolys->InsertCellPoint(id0);
			newPolys->InsertCellPoint(id1);
			newPolys->InsertCellPoint(id2);
		}
		else
		{
			NnonManifolds++;
		}
	}

	output->SetPoints(pts);
	output->SetPolys(newPolys);
	newPolys->Delete();
	pts->Delete();
	neighbors->Delete();

	std::cout << "Found and removed " << NnonManifolds << " triangles with non-manifold edges" << std::endl;

	return 1;
}

//--------------------------------------------------------------------------
void vtkRemoveNonManifoldCells::PrintSelf(ostream& os, vtkIndent indent) 
{
	this->Superclass::PrintSelf(os,indent);
}

//--------------------------------------------------------------------------
vtkMTimeType vtkRemoveNonManifoldCells::GetMTime()
{
	unsigned long mTime=this->vtkObject::GetMTime();
	return mTime;
}
