#include "vtkPolyDataMarkCoveredCells.h"

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

vtkStandardNewMacro(vtkPolyDataMarkCoveredCells);

vtkPolyDataMarkCoveredCells::vtkPolyDataMarkCoveredCells()
{
	MarkRadius = 1.0;
	Source = NULL;
}

//--------------------------------------------------------------------------
vtkPolyDataMarkCoveredCells::~vtkPolyDataMarkCoveredCells()
{
}


//--------------------------------------------------------------------------
int vtkPolyDataMarkCoveredCells::RequestData(
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

	vtkDebugMacro(<<"Beginning vtkPolyDataMarkCoveredCells");
	std::cout << "Beginning vtkPolyDataMarkCoveredCells" << std::endl;
	if ( (numPts<1) || (inPts == NULL ) )
	{
		vtkDebugMacro(<<"No data to Operate On!");
		return 1;
	}
	if (Source == NULL)
	{
		vtkDebugMacro(<<"Source not set!");
		return 1;
	}

	output->DeepCopy(input);

	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);
	scalars->SetNumberOfValues(input->GetNumberOfPoints());

	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		scalars->SetValue(i, 0);
	}

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	double markval = 1;
	for (int i = 0; i < Source->GetNumberOfPoints(); i++)
	{
		double p[3];
		Source->GetPoint(i, p);

		vtkIdList *idlist = vtkIdList::New();
		locator->FindPointsWithinRadius(MarkRadius, p, idlist);

		for (int j = 0; j < idlist->GetNumberOfIds(); j++)
		{
			int id = idlist->GetId(j);

			scalars->SetValue(id, markval);
		}
		idlist->Delete();
	}

	locator->Delete();

	output->GetPointData()->SetScalars(scalars);

	scalars->Delete();

	return 1;
}

//--------------------------------------------------------------------------
void vtkPolyDataMarkCoveredCells::PrintSelf(ostream& os, vtkIndent indent) 
{
	this->Superclass::PrintSelf(os,indent);
}

//--------------------------------------------------------------------------
vtkMTimeType vtkPolyDataMarkCoveredCells::GetMTime()
{
	unsigned long mTime=this->vtkObject::GetMTime();
	return mTime;
}
