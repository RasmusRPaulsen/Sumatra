#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkPolyDataProjection.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include <vtkTriangle.h>
#include <iostream>


vtkStandardNewMacro(vtkPolyDataProjection);

//----------------------------------------------------------------------------
vtkPolyDataProjection::vtkPolyDataProjection()
{
	CreateDisplacementField = false;
	m_DisplacementField = vtkPolyData::New();
	this->SetNumberOfInputPorts(2);
}
//----------------------------------------------------------------------------
vtkPolyDataProjection::~vtkPolyDataProjection()
{
	if (m_DisplacementField)
		m_DisplacementField->Delete();
}

void vtkPolyDataProjection::SetTargetData(vtkPolyData *target)
{
  this->SetInputData(1, target);
}

vtkPolyData *vtkPolyDataProjection::GetTarget()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkPolyData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
}

int vtkPolyDataProjection::FillInputPortInformation(int port,
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


//----------------------------------------------------------------------------
int vtkPolyDataProjection::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
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
	
	output->CopyStructure(input);
	output->GetPointData()->PassData(input->GetPointData());
	output->GetCellData()->PassData(input->GetCellData());

	vtkCellLocator *locator = vtkCellLocator::New();

	locator->SetDataSet(target);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();

	vtkPoints *pts = vtkPoints::New();
	pts->SetNumberOfPoints(input->GetNumberOfPoints());

	vtkCellArray *verts = vtkCellArray::New();


	for (int k = 0; k < input->GetPoints()->GetNumberOfPoints(); k++)
	{
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
	
		output->GetPoints()->SetPoint(k, cp);

		double disp[3];
		disp[0] = cp[0] - P[0];
		disp[1] = cp[1] - P[1];
		disp[2] = cp[2] - P[2];
		
		pts->SetPoint(k, disp);

		verts->InsertNextCell(1);
		verts->InsertCellPoint(k);
	}

	m_DisplacementField->SetPoints(pts);
	pts->Delete();
	m_DisplacementField->SetVerts(verts);
	verts->Delete();

	vtkDebugMacro(<< "<--Ending vtkPolyDataProjection");
	return 1;
}


void vtkPolyDataProjection::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

vtkPolyData * vtkPolyDataProjection::GetDisplacementField()
{
	return m_DisplacementField;
}