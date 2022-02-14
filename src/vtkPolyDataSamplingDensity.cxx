#include "vtkPolyDataSamplingDensity.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
//#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
//#include "vtkPolygon.h"
//#include "vtkTriangleStrip.h"
//#include "vtkPriorityQueue.h"

#include <vtkDoubleArray.h>
#include <vtkPointLocator.h>

vtkStandardNewMacro(vtkPolyDataSamplingDensity);


vtkPolyDataSamplingDensity::vtkPolyDataSamplingDensity()
{
	SearchRadius = 1.0;
}

void vtkPolyDataSamplingDensity::ComputeLocalDensity(vtkPolyData *input, vtkPolyData *output)
{
	std::cout << "Computing local density" << std::endl;

	double V = 4 * vtkMath::Pi() * SearchRadius * SearchRadius * SearchRadius;

	vtkPointLocator *Locator = vtkPointLocator::New();
	Locator->SetDataSet(input);
	Locator->SetNumberOfPointsPerBucket(1);
	Locator->BuildLocator();

	vtkDoubleArray *outscalars = vtkDoubleArray::New();
	outscalars->SetNumberOfComponents(1);

	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double p[3];
		input->GetPoint(i, p);

		vtkIdList *idlist = vtkIdList::New();
		Locator->FindPointsWithinRadius(SearchRadius, p, idlist);

		int N = idlist->GetNumberOfIds();

		double density = (double)N / V;

		outscalars->InsertNextTuple(&density);

		idlist->Delete();
	}

	Locator->Delete();

	output->GetPointData()->SetScalars(outscalars);

	outscalars->Delete();
}

// Generate normals for polygon meshes
int vtkPolyDataSamplingDensity::RequestData(
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

	output->DeepCopy(input);

	ComputeLocalDensity(input, output);

	return 1;
}

void vtkPolyDataSamplingDensity::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
