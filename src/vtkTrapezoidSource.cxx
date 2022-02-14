#include "vtkTrapezoidSource.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include <vtkMath.h>

#include <math.h>

vtkStandardNewMacro(vtkTrapezoidSource);

vtkTrapezoidSource::vtkTrapezoidSource(double xL, double yL, double zL, double A, double B, int Cap)
{
  this->XLength = fabs(xL);
  this->YLength = fabs(yL);
  this->ZLength = fabs(zL);
  this->Alpha   = A;
  this->Beta    = B;
  this->Capping = Cap;

  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;

  this->SetNumberOfInputPorts(0);
}

int vtkTrapezoidSource::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  double x[3], n[3];
  int numPolys=6, numPts=24;
  vtkIdType pts[4];
  vtkPoints *newPoints; 
  vtkFloatArray *newNormals;
  vtkCellArray *newPolys;

  if (!Capping)
  {
	  numPolys = 4;
	  numPts   = 16;
  }

  newPoints = vtkPoints::New();
  newPoints->Allocate(numPts);
  newNormals = vtkFloatArray::New();
  newNormals->SetNumberOfComponents(3);
  newNormals->Allocate(numPts);
  newNormals->SetName("Normals");
  
  newPolys = vtkCellArray::New();
  newPolys->Allocate(newPolys->EstimateSize(numPolys,4));

	// Create bottom poly. (ok)
	n[0] = 0; n[1] = -1; n[2] = 0;

	x[0] = Center[0] - XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] - XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);
	
	pts[0] = 0; pts[1] = 1; pts[2] = 2; pts[3] = 3; 
	newPolys->InsertNextCell(4,pts);

	// Top poly (ok)
	double dA = tan(vtkMath::RadiansFromDegrees(90.0 - Alpha)) * YLength;
	double dB = tan(vtkMath::RadiansFromDegrees(90.0 - Beta)) * YLength;
	n[0] = 0; n[1] =  1; n[2] = 0;

	x[0] = Center[0] - XLength / 2.0 + dA;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] - XLength / 2.0 + dA;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0 - dB;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0 - dB;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);
	
	pts[0] = 4; pts[1] = 5; pts[2] = 6; pts[3] = 7; 
	newPolys->InsertNextCell(4,pts);

	// Side poly
	n[0] = -sin(vtkMath::RadiansFromDegrees(Alpha)); n[1] =  cos(vtkMath::RadiansFromDegrees(Alpha)); n[2] = 0;
	
	x[0] = Center[0] - XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] - XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] - XLength / 2.0 + dA;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] - XLength / 2.0 + dA;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	pts[0] = 8; pts[1] = 9; pts[2] = 10; pts[3] = 11; 
	newPolys->InsertNextCell(4,pts);

	// Side poly 2
	n[0] = sin(vtkMath::RadiansFromDegrees(Beta)); n[1] =  cos(vtkMath::RadiansFromDegrees(Beta)); n[2] = 0;
	
	x[0] = Center[0] + XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0;
	x[1] = Center[1] - YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0 - dB;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] - ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	x[0] = Center[0] + XLength / 2.0 - dB;
	x[1] = Center[1] + YLength / 2.0;
	x[2] = Center[2] + ZLength / 2.0;
	newPoints->InsertNextPoint(x);
	newNormals->InsertNextTuple(n);

	pts[0] = 12; pts[1] = 13; pts[2] = 14; pts[3] = 15; 
	newPolys->InsertNextCell(4,pts);

	if (Capping)
	{
		// End poly (ok)
		n[0] = 0; n[1] =  0; n[2] = 1; 
		x[0] = Center[0] - XLength / 2.0;
		x[1] = Center[1] - YLength / 2.0;
		x[2] = Center[2] + ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		x[0] = Center[0] + XLength / 2.0;
		x[1] = Center[1] - YLength / 2.0;
		x[2] = Center[2] + ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		x[0] = Center[0] + XLength / 2.0 - dB;
		x[1] = Center[1] + YLength / 2.0;
		x[2] = Center[2] + ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		x[0] = Center[0] - XLength / 2.0 + dA;
		x[1] = Center[1] + YLength / 2.0;
		x[2] = Center[2] + ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		pts[0] = 16; pts[1] = 17; pts[2] = 18; pts[3] = 19; 
		newPolys->InsertNextCell(4,pts);

		// End poly 2 (ok)
		n[0] = 0; n[1] =  0; n[2] = -1; 
		x[0] = Center[0] + XLength / 2.0;
		x[1] = Center[1] - YLength / 2.0;
		x[2] = Center[2] - ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		x[0] = Center[0] - XLength / 2.0;
		x[1] = Center[1] - YLength / 2.0;
		x[2] = Center[2] - ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		x[0] = Center[0] - XLength / 2.0 + dA;
		x[1] = Center[1] + YLength / 2.0;
		x[2] = Center[2] - ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		x[0] = Center[0] + XLength / 2.0 - dB;
		x[1] = Center[1] + YLength / 2.0;
		x[2] = Center[2] - ZLength / 2.0;
		newPoints->InsertNextPoint(x);
		newNormals->InsertNextTuple(n);

		pts[0] = 20; pts[1] = 21; pts[2] = 22; pts[3] = 23; 
		newPolys->InsertNextCell(4,pts);
	}

// Update ourselves and release memory
//
  output->SetPoints(newPoints);
  newPoints->Delete();

  output->GetPointData()->SetNormals(newNormals);
  newNormals->Delete();

  newPolys->Squeeze(); // since we've estimated size; reclaim some space
  output->SetPolys(newPolys);
  newPolys->Delete();

  return 1;
}

void vtkTrapezoidSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "X Length: " << this->XLength << "\n";
  os << indent << "Y Length: " << this->YLength << "\n";
  os << indent << "Z Length: " << this->ZLength << "\n";
  os << indent << "Alpha: " << this->Alpha << "\n";
  os << indent << "Beta: " << this->Beta << "\n";
  os << indent << "Center: (" << this->Center[0] << ", " 
               << this->Center[1] << ", " << this->Center[2] << ")\n";
}
