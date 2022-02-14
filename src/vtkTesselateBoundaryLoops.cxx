/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTesselateBoundaryLoops.cxx
  Language:  C++
  Version:   $Id: vtkTesselateBoundaryLoops.cxx 1.1 2004/04/16 11:07:16 Goodwin Exp Goodwin $

  Copyright (c) 2004 Goodwin Lawlor
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTesselateBoundaryLoops.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"


vtkStandardNewMacro(vtkTesselateBoundaryLoops);

vtkTesselateBoundaryLoops::vtkTesselateBoundaryLoops()
{
  this->Tolerance = 0.0;
  this->AppendTesselationToInput = 0;
  this->ReverseSense = 0;
  this->MaximumLength = 1000;
  this->NumberOfLoops = 0;
  this->NonManifoldEdges = 0;
  
  // setup the pipeline
  this->FeatureEdges = vtkFeatureEdges::New();
    this->FeatureEdges->BoundaryEdgesOn();
    this->FeatureEdges->FeatureEdgesOff();
  this->Clean  = vtkCleanPolyData::New();
    this->Clean->SetInputConnection(this->FeatureEdges->GetOutputPort());
  this->Stripper = vtkStripper::New();
    this->Stripper->SetInputConnection(Clean->GetOutputPort());
  this->tmpPolyData = vtkPolyData::New();
  this->Reverse = vtkReverseSense::New();
    this->Reverse->SetInputData(tmpPolyData);
  this->Tesselator = vtkTriangleFilter::New();
    this->Tesselator->SetInputConnection(this->Reverse->GetOutputPort());
  this->Append = vtkAppendPolyData::New();
    this->Append->AddInputConnection(this->Tesselator->GetOutputPort());
  this->Clean2 = vtkCleanPolyData::New();
      this->Clean2->SetInputConnection(this->Append->GetOutputPort());
}

vtkTesselateBoundaryLoops::~vtkTesselateBoundaryLoops()
{
  this->FeatureEdges->Delete();
  this->Clean->Delete();
  this->Stripper->Delete();
  this->tmpPolyData->Delete();
  this->Reverse->Delete();
  this->Tesselator->Delete();
  this->Append->Delete();
  this->Clean2->Delete();
}


int vtkTesselateBoundaryLoops::RequestData(
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
  
  
  this->FeatureEdges->SetInputData(input);
  // sometimes a boundary edge is non-manifold.
  this->FeatureEdges->SetNonManifoldEdges(this->NonManifoldEdges);    

  this->Clean->SetTolerance(this->Tolerance);
  
  this->Stripper->SetMaximumLength(this->MaximumLength);
  this->Stripper->Update();
  this->NumberOfLoops = this->Stripper->GetOutput()->GetNumberOfLines();
  
  this->tmpPolyData->SetPolys(this->Stripper->GetOutput()->GetLines());
  this->tmpPolyData->SetPoints(this->Stripper->GetOutput()->GetPoints());
  
  this->Reverse->SetReverseCells(this->ReverseSense);
  
  vtkPolyData *last;
  
  if (this->AppendTesselationToInput)
    {
    this->Append->AddInputData(input);
    this->Clean2->Update();
    last = this->Clean2->GetOutput();
    }
  else
    {
    this->Tesselator->Update();
    last = this->Tesselator->GetOutput();
    }
  
  output->SetPoints(last->GetPoints());
  output->SetPolys(last->GetPolys());
  
  vtkDebugMacro(<< "Finished tesselation");

  return 1;
  
}

void vtkTesselateBoundaryLoops::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}