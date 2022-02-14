#include "vtkRemoveUnusedPolyDataPoints.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <iostream>
#include <set>
#include <map>

vtkStandardNewMacro(vtkRemoveUnusedPolyDataPoints);

//----------------------------------------------------------------------------
vtkRemoveUnusedPolyDataPoints::vtkRemoveUnusedPolyDataPoints()
{

}

vtkRemoveUnusedPolyDataPoints::~vtkRemoveUnusedPolyDataPoints()
{
}

int vtkRemoveUnusedPolyDataPoints::RequestData(
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
  

	vtkPointData *pd = input->GetPointData();
	vtkCellData *cd = input->GetCellData();
	vtkPointData *outputPD = output->GetPointData();
	vtkCellData *outputCD = output->GetCellData();
	vtkPoints *inPts=input->GetPoints();
	vtkIdType numPts, i;
	vtkCellArray *inVerts=NULL, *inLines=NULL, *inPolys=NULL, *inStrips=NULL;
	vtkCellArray *newVerts=NULL, *newLines=NULL, *newPolys=NULL, *newStrips=NULL;
	

	numPts = input->GetNumberOfPoints();
	
//	output->SetPoints(inPts);
//	outputPD->PassData(pd);
	
	std::set<vtkIdType> usedIDs;
	
	// Now loop over all cells to see whether they have points that are in the id list
	// Copy if they have not. Note: there is an awful hack here, that
	// can result in bugs. The cellId is assumed to be arranged starting
	// with the verts, then lines, then polys, then strips.
	//
//	int numIn;
	vtkIdType npts;
	const vtkIdType *pts;
	if ( input->GetNumberOfVerts() )
	{
		inVerts = input->GetVerts();
		newVerts = vtkCellArray::New();
		newVerts->Allocate(inVerts->GetSize());
	}
	if ( input->GetNumberOfLines() )
	{
		inLines = input->GetLines();
		newLines = vtkCellArray::New();
		newLines->Allocate(inLines->GetSize());
	}
	if ( input->GetNumberOfPolys() )
	{
		inPolys = input->GetPolys();
		newPolys = vtkCellArray::New();
		newPolys->Allocate(inPolys->GetSize());
	}
	if ( input->GetNumberOfStrips() )
	{
		inStrips = input->GetStrips();
		newStrips = vtkCellArray::New();
		newStrips->Allocate(inStrips->GetSize());
	}
	
	// verts
	if ( newVerts && !this->GetAbortExecute() )
	{
		for (inVerts->InitTraversal(); inVerts->GetNextCell(npts,pts); )
		{
			for (i=0; i<npts; i++)
			{
				usedIDs.insert(pts[i]);
			}
		}
	}
	this->UpdateProgress (0.1);
	
	// lines
	if ( newLines && !this->GetAbortExecute() )
	{
		for (inLines->InitTraversal(); inLines->GetNextCell(npts,pts); )
		{
			for (i=0; i<npts; i++)
			{
				usedIDs.insert(pts[i]);
			}
		}
	}
	this->UpdateProgress (0.2);
	
	// polys
	if ( newPolys && !this->GetAbortExecute() )
	{
		for (inPolys->InitTraversal(); inPolys->GetNextCell(npts,pts); )
		{
			for (i=0; i<npts; i++)
			{
				usedIDs.insert(pts[i]);
			}
		}
	}
	this->UpdateProgress (0.3);
	
	// strips
	if ( newStrips && !this->GetAbortExecute() )
	{
		for (inStrips->InitTraversal(); inStrips->GetNextCell(npts,pts); )
		{
			for (i=0; i<npts; i++)
			{
				usedIDs.insert(pts[i]);
			}
		}
	}
	this->UpdateProgress (0.4);
	
	std::cout << "UsedIDs size: " << usedIDs.size() << " out of " << input->GetNumberOfPoints() << " input points" << std::endl;

	// Create a mapping from the old used id values to new ids
	std::map<vtkIdType, vtkIdType> idLUT;

	vtkIdType newId = 0;
	std::set<vtkIdType>::iterator it;
	for (it = usedIDs.begin(); it != usedIDs.end(); it++)
	{
		idLUT.insert(std::pair<vtkIdType, vtkIdType>(*it, newId));
//		std::cout << *it << ", ";
		newId++;
	}

	std::map<vtkIdType, vtkIdType>::iterator mapIt;
//	mapIt = idLUT.find(3);
//	std::cout << "test value " << (*mapIt).second << std::endl;

	size_t nNewPoints = usedIDs.size();
	vtkPoints *newPts = vtkPoints::New();
	 newPts->SetNumberOfPoints(nNewPoints);

	// Create new point array
	for (mapIt = idLUT.begin(); mapIt != idLUT.end(); mapIt++)
	{
		vtkIdType oldId = (*mapIt).first;
		vtkIdType newId = (*mapIt).second;
		newPts->SetPoint(newId, inPts->GetPoint(oldId));
	}

	vtkDoubleArray *newscals = NULL;

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetScalars());
	if (scalars)
	{
		newscals = vtkDoubleArray::New();
		newscals->SetNumberOfComponents(1);
		newscals->SetNumberOfValues(nNewPoints);

		for (mapIt = idLUT.begin(); mapIt != idLUT.end(); mapIt++)
		{
			vtkIdType oldId = (*mapIt).first;
			vtkIdType newId = (*mapIt).second;
			double val = scalars->GetValue(oldId);
			newscals->SetValue(newId, val);
		}
	}

	vtkDoubleArray *newnorms = NULL;

	vtkDataArray *normals = input->GetPointData()->GetNormals();
	if (normals)
	{
		newnorms = vtkDoubleArray::New();
		newnorms->SetNumberOfComponents(3);
		newnorms->SetName("Normals");
		newnorms->SetNumberOfTuples(nNewPoints);

		for (mapIt = idLUT.begin(); mapIt != idLUT.end(); mapIt++)
		{
			vtkIdType oldId = (*mapIt).first;
			vtkIdType newId = (*mapIt).second;
			double n[3];
			normals->GetTuple(oldId, n);
			newnorms->SetTuple(newId, n);
		}
	}


	this->UpdateProgress (0.5);

	// update new verts
	if ( newVerts && !this->GetAbortExecute() )
	{
		for (inVerts->InitTraversal(); inVerts->GetNextCell(npts,pts); )
		{
	        vtkIdType newId = newVerts->InsertNextCell(npts);

			// Translate point ids
			for (i=0; i<npts; i++)
			{
				mapIt = idLUT.find(pts[i]);
				vtkIdType newId = (*mapIt).second;
				newVerts->InsertCellPoint(newId);
			}
		}
	}

	// update new lines
	if ( newLines && !this->GetAbortExecute() )
	{
		for (inLines->InitTraversal(); inLines->GetNextCell(npts,pts); )
		{
	        vtkIdType newId = newLines->InsertNextCell(npts);

			// Translate point ids
			for (i=0; i<npts; i++)
			{
				mapIt = idLUT.find(pts[i]);
				vtkIdType newId = (*mapIt).second;
				newLines->InsertCellPoint(newId);
			}
		}
	}

	this->UpdateProgress (0.6);
	// update new polys
	if ( newPolys && !this->GetAbortExecute() )
	{
		for (inPolys->InitTraversal(); inPolys->GetNextCell(npts,pts); )
		{
	        vtkIdType newId = newPolys->InsertNextCell(npts);

			// Translate point ids
			for (i=0; i<npts; i++)
			{
				mapIt = idLUT.find(pts[i]);
				vtkIdType newId = (*mapIt).second;
				newPolys->InsertCellPoint(newId);
			}
		}
	}

	// update new strips
	if ( newStrips && !this->GetAbortExecute() )
	{
		for (inStrips->InitTraversal(); inStrips->GetNextCell(npts,pts); )
		{
	        vtkIdType newId = newStrips->InsertNextCell(npts);

			// Translate point ids
			for (i=0; i<npts; i++)
			{
				mapIt = idLUT.find(pts[i]);
				vtkIdType newId = (*mapIt).second;
				newStrips->InsertCellPoint(newId);
			}
		}
	}

	this->UpdateProgress (0.7);

	// Update ourselves and release memory
	//
	if (newPts)
	{
		output->SetPoints(newPts);
		newPts->Delete();
	}	
	if ( newVerts )
	{
		output->SetVerts(newVerts);
		newVerts->Delete();
	}
	if ( newLines )
	{
		output->SetLines(newLines);
		newLines->Delete();
	}
	if ( newPolys )
	{
		output->SetPolys(newPolys);
		newPolys->Delete();
	}
	if ( newStrips )
	{
		output->SetStrips(newStrips);
		newStrips->Delete();
	}
	if (newscals)
	{
		output->GetPointData()->SetScalars(newscals);
		newscals->Delete();
	}
	if (newnorms)
	{
		output->GetPointData()->SetNormals(newnorms);
		newnorms->Delete();
	}
	this->UpdateProgress (1.0);
	return 1;
}

void vtkRemoveUnusedPolyDataPoints::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
