#include "vtkPolyDataResamplePoints.h"

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

vtkStandardNewMacro(vtkPolyDataResamplePoints);


vtkPolyDataResamplePoints::vtkPolyDataResamplePoints()
{
	NumberOfOutputPoints = 1000000;
}

// Generate normals for polygon meshes
int vtkPolyDataResamplePoints::RequestData(
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

	int Nin = input->GetNumberOfPoints();
	if (Nin < 1)
		return 0;

	if (NumberOfOutputPoints > Nin)
	{
		std::cout << "Input points " << Nin << " you asked for " << NumberOfOutputPoints << " input points are copied" << std::endl;
		output->DeepCopy(input);
		return 1;
	}


	double step = Nin / (double)NumberOfOutputPoints;

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *Normals = vtkDoubleArray::New();
	Normals->SetNumberOfComponents(3);
	Normals->SetName("Normals");
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	vtkDataArray *innorms = input->GetPointData()->GetNormals(); 
	vtkDoubleArray *inscalars = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetScalars());


	double pos = 0;
	for (int i = 0; i < NumberOfOutputPoints; i++)
	{
		double p[3];
		int ipos = (int)pos;
		input->GetPoint(ipos, p);
		vtkIdType id = pts->InsertNextPoint(p);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);

		if (innorms)
		{
			double n[3];
			innorms->GetTuple(ipos, n);
			Normals->InsertNextTuple(n);
		}
		if (inscalars)
		{
			double val = inscalars->GetValue(ipos);
			scalars->InsertNextValue(val);
		}
		pos += step;
		if (pos >= Nin)
			pos = pos - Nin;
	}

	output->SetPoints(pts);
	output->SetVerts(verts);
	if (innorms)
	{
		output->GetPointData()->SetNormals(Normals);
	}
	if (inscalars)
	{
		output->GetPointData()->SetScalars(scalars);
	}

	Normals->Delete();
	pts->Delete();
	verts->Delete();
	scalars->Delete();

	return 1;
}

void vtkPolyDataResamplePoints::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
