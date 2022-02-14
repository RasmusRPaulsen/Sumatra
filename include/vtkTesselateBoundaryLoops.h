/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTesselateBoundaryLoops.h
  Language:  C++
  Version:   $Id: vtkTesselateBoundaryLoops.h 1.2 2004/07/09 16:02:12 Goodwin Exp Goodwin $

  Copyright (c) 2004 Goodwin Lawlor
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTesselateBoundaryLoops - tesselate boundary loops, aka, fill in surface holes.
// .SECTION Description
// This filter will convert boundary edges that form a  closed loop into tesselated polygons. It will not
// smoothly interpolate over a surface hole.



#ifndef __vtkTesselateBoundaryLoops_h
#define __vtkTesselateBoundaryLoops_h

#include "vtkPolyDataAlgorithm.h"
//#include "vtkBioengConfigure.h" // Include configuration header.
#include "vtkFeatureEdges.h"
#include "vtkCleanPolyData.h"
#include "vtkStripper.h"
#include "vtkPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkReverseSense.h"

class  vtkTesselateBoundaryLoops : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkTesselateBoundaryLoops,vtkPolyDataAlgorithm);
  static vtkTesselateBoundaryLoops *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description:
  // Set the tolerance for point merging of line verts... if line verts aren't exactly coincident
  // they wont be joined up into a loop. With a tolerance they can be incouraged to be coincident.
  vtkSetMacro(Tolerance, float);
  vtkGetMacro(Tolerance, float);
  
  // Description:
  // Set the maximum number of lines in a loop. Default is 1000;
  vtkSetMacro(MaximumLength, int);
  vtkGetMacro(MaximumLength, int);
  
  //Description:
  // Append the tesselated loops to the input dataset. Default is Off.
  vtkSetMacro(AppendTesselationToInput, int);
  vtkGetMacro(AppendTesselationToInput, int);
  vtkBooleanMacro(AppendTesselationToInput,int);
  
  // Description:
  // Reverse the sense (front/back face) of the loops being tesselated
  vtkSetMacro(ReverseSense, int);
  vtkGetMacro(ReverseSense, int);
  vtkBooleanMacro(ReverseSense, int);
  
  // Description:
  // Sometimes boundary edges are non-manifold and need to be included
  vtkSetMacro(NonManifoldEdges, int);
  vtkGetMacro(NonManifoldEdges, int);
  vtkBooleanMacro(NonManifoldEdges, int);
  
  // Description:
  // Get the number of boundary loops 
  vtkGetMacro(NumberOfLoops, int);
  
protected:
  vtkTesselateBoundaryLoops();
  ~vtkTesselateBoundaryLoops();

  vtkFeatureEdges *FeatureEdges;
  vtkCleanPolyData *Clean;
  vtkStripper *Stripper;
  vtkPolyData *tmpPolyData;
  vtkReverseSense *Reverse;
  vtkTriangleFilter *Tesselator;
  vtkAppendPolyData *Append;
  vtkCleanPolyData *Clean2;
  
  float Tolerance;
  int AppendTesselationToInput;
  int ReverseSense;
  int MaximumLength;
  int NumberOfLoops;
  int NonManifoldEdges;
  
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);



private:
  vtkTesselateBoundaryLoops(const vtkTesselateBoundaryLoops&);  // Not implemented.
  void operator=(const vtkTesselateBoundaryLoops&);  // Not implemented.

};

#endif

