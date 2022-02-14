/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPrincipalAxisTransform.h,v $
  Language:  C++
  Date:      $Date: 2002/04/18 13:52:49 $
  Version:   $Revision: 1.5 $
  
  Made by Rasmus Paulsen
  email:  rrp@imm.dtu.dk
  web:    www.imm.dtu.dk/~rrp/VTK

  This class is not mature enough to enter the official VTK release.
=========================================================================*/
// .NAME vtkPrincipalAxisTransform - Create a transform that scales, translates and rotates an object to lie on the input datas principal axis.
// .SECTION Description
//
// vtkPrincipalAxisTransform generates a transform that scales, translates and rotates
// a given object to the principal axis of the input data.
// If the input data is a Gaussian distributed point cloud and the transformation
// is applied to a unit sphere, the transformed sphere will show the covariance structure
// of the point cloud.
// It is possible to specify if the transform shall include 
// rotation, translation, scaling or all of them.

#ifndef __vtkPrincipalAxisTransform_h
#define __vtkPrincipalAxisTransform_h

#include "vtkTransform.h"

class vtkPoints;

class vtkPrincipalAxisTransform : public vtkTransform
{
public:
  static vtkPrincipalAxisTransform *New();
  vtkTypeMacro(vtkPrincipalAxisTransform,vtkLinearTransform);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the sourcedata sets.
  void SetSource(vtkPoints *source);
  vtkGetObjectMacro(Source, vtkPoints);

  // Description: 
  // Use scale
  vtkSetMacro(DoScale, int);
  vtkGetMacro(DoScale, int);
  vtkBooleanMacro(DoScale, int);

  // Description: 
  // Use Translation
  vtkSetMacro(DoTranslate, int);
  vtkGetMacro(DoTranslate, int);
  vtkBooleanMacro(DoTranslate, int);

  // Description: 
  // Use Rotation
  vtkSetMacro(DoRotate, int);
  vtkGetMacro(DoRotate, int);
  vtkBooleanMacro(DoRotate, int);

  // Description:
  // Make another transform of the same type.
  vtkAbstractTransform *MakeTransform();

protected:

  // Description:
  // Release source
  void ReleaseSource(void);

  // Description:
  // Get the MTime of this object
  vtkMTimeType GetMTime();

  vtkPrincipalAxisTransform();
  ~vtkPrincipalAxisTransform();
  vtkPrincipalAxisTransform(const vtkPrincipalAxisTransform&) {};
  void operator=(const vtkPrincipalAxisTransform&) {};

  void InternalUpdate();

  // Description:
  // This method does no type checking, use DeepCopy instead.
  void InternalDeepCopy(vtkAbstractTransform *transform);

  vtkPoints* Source;

  int DoScale;

  int DoTranslate;

  int DoRotate;
};

#endif
