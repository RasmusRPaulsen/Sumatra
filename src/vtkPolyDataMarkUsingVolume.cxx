#include "vtkPolyDataMarkUsingVolume.h"


#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include <vtkImplicitVolume.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>


vtkStandardNewMacro(vtkPolyDataMarkUsingVolume);

vtkPolyDataMarkUsingVolume::vtkPolyDataMarkUsingVolume()
{
	Volume = NULL;
}

//--------------------------------------------------------------------------
vtkPolyDataMarkUsingVolume::~vtkPolyDataMarkUsingVolume()
{
}


//--------------------------------------------------------------------------
int vtkPolyDataMarkUsingVolume::RequestData(
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

	vtkDebugMacro(<<"Beginning vtkPolyDataMarkUsingVolume");
	std::cout << "Beginning vtkPolyDataMarkUsingVolume" << std::endl;
	if ( (numPts<1) || (inPts == NULL ) )
	{
		vtkDebugMacro(<<"No data to Operate On!");
		return 1;
	}
	if (Volume == NULL)
	{
		vtkDebugMacro(<<"No volume set!");
		return 1;
	}

	output->DeepCopy(input);

	vtkImplicitVolume *impvol = vtkImplicitVolume::New();
	impvol->SetVolume(Volume);

	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);
	scalars->SetNumberOfValues(input->GetNumberOfPoints());

	for (int i = 0; i < output->GetNumberOfPoints(); i++)
	{
		double p[3];
		output->GetPoint(i, p);

		double val = impvol->EvaluateFunction(p);

		if (val < 0)
		{
			val = 0;
//			std::cout << "Trouble at " << i << " (" << p[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;
		}
		if (val > 1.0)
		{
			val = 1.0;
		}

		scalars->SetValue(i, val);
	}
	output->GetPointData()->SetScalars(scalars);

	scalars->Delete();
	impvol->Delete();

	return 1;
}


//--------------------------------------------------------------------------
void vtkPolyDataMarkUsingVolume::PrintSelf(ostream& os, vtkIndent indent) 
{
	this->Superclass::PrintSelf(os,indent);
}

//--------------------------------------------------------------------------
vtkMTimeType vtkPolyDataMarkUsingVolume::GetMTime()
{
	unsigned long mTime=this->vtkObject::GetMTime();
	return mTime;
}
