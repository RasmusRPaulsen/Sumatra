#include "vtkCropPolyData.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkMergePoints.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include <vtkPointLocator.h>
#include <vtkDoubleArray.h>

#include <deque>

vtkStandardNewMacro(vtkCropPolyData);

vtkCropPolyData::vtkCropPolyData()
{
	Tolerance = 1e-6;
}

//--------------------------------------------------------------------------
vtkCropPolyData::~vtkCropPolyData()
{
}


//--------------------------------------------------------------------------
int vtkCropPolyData::RequestData(
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

	vtkDebugMacro(<<"Beginning vtkCropPolyData");
	std::cout << "Beginning vtkCropPolyData" << std::endl;
	if ( (numPts<1) || (inPts == NULL ) )
	{
		vtkDebugMacro(<<"No data to Operate On!");
		return 1;
	}

	output->DeepCopy(input);
//	PointStatus.resize(input->GetNumberOfPoints(), 0);

	LocateEdgePoints(input, output);
	MarkConnectedPoints(output);

	return 1;
}
void vtkCropPolyData::LocateEdgePoints(vtkPolyData *input, vtkPolyData *output)
{
//	double tolerance = 1e-6;
//	double tolerance = 5;

	//double bounds[6];
	//input->GetBounds(bounds);

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetScalars());
	if (!scalars)
	{
		std::cerr << "No scalars in input" << std::endl;
		return;
	}
	vtkDoubleArray *outscalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
	if (!outscalars)
	{
		std::cerr << "No scalars in output" << std::endl;
		return;
	}


	int edgep = 0;
	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double p[3];
		input->GetPoint(i, p);

		if (scalars->GetValue(i) == 0 && 
			(abs(p[0]-Boundary[0]) < Tolerance || abs(p[0]-Boundary[1]) < Tolerance ||
			abs(p[1]-Boundary[2]) < Tolerance || abs(p[1]-Boundary[3]) < Tolerance ||
			abs(p[2]-Boundary[4]) < Tolerance || abs(p[2]-Boundary[5]) < Tolerance ))
		{
//			PointStatus[i] = 1;
			outscalars->SetValue(i, 2);
			edgep++;
		}
	}
	outscalars->Modified();

	std::cout << "Found " << edgep << " edgepoints" << std::endl;
}


void vtkCropPolyData::MarkConnectedPoints( vtkPolyData *output )
{
	std::deque<vtkIdType> visitlist;

	vtkDoubleArray *outscalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
	if (!outscalars)
	{
		std::cerr << "No scalars in output" << std::endl;
		return;
	}

	output->BuildLinks();

	// Push seeds on visitlist
	for (int i = 0; i < output->GetNumberOfPoints(); i++)
	{
		if (outscalars->GetValue(i) == 2)
		{
			visitlist.push_back(i);
		}
	}
	vtkIdList *verts  = vtkIdList::New();
	vtkIdList *cellids = vtkIdList::New();

	while (!visitlist.empty())
	{
		vtkIdType curId = visitlist.front(); 
		visitlist.pop_front();

		outscalars->SetValue(curId, 2);

		// Get the cells that point with curID belongs to:
		output->GetPointCells(curId, cellids);

		for (int nc = 0; nc < cellids->GetNumberOfIds(); nc++)
		{
			// Cell ID
			vtkIdType cid = cellids->GetId(nc);

			// Get the points of the neighbour cell
			// it is assumed that all points in a cell is connected...
			output->GetCellPoints(cid, verts);

			// Iterate over all neigbhour points
			for (int nv = 0; nv < verts->GetNumberOfIds(); nv++)
			{
				vtkIdType tid = verts->GetId(nv);

				if (outscalars->GetValue(tid) == 0)
				{
					outscalars->SetValue(tid, 2);

					visitlist.push_back(tid);
				}
			}
		}
	}
	verts->Delete();
	cellids->Delete();
	outscalars->Modified();
}


//--------------------------------------------------------------------------
void vtkCropPolyData::PrintSelf(ostream& os, vtkIndent indent) 
{
	this->Superclass::PrintSelf(os,indent);
}

//--------------------------------------------------------------------------
vtkMTimeType vtkCropPolyData::GetMTime()
{
	unsigned long mTime=this->vtkObject::GetMTime();
	return mTime;
}

void vtkCropPolyData::SetBoundary( double boundary[6] )
{
	for (int i = 0; i < 6; i++)
	{
		Boundary[i] = boundary[i];
	}
}