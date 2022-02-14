/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkLandmarkTransformDirectQuart.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkLandmarkTransformDirectQuart - a linear transform specified by two corresponding point sets
// .SECTION Description
// A vtkLandmarkTransformDirectQuart is defined by two sets of landmarks, the 
// transform computed gives the best fit mapping one onto the other, in a 
// least squares sense. The indices are taken to correspond, so point 1 
// in the first set will get mapped close to point 1 in the second set, 
// etc. Call SetSourceLandmarks and SetTargetLandmarks to specify the two
// sets of landmarks, ensure they have the same number of points.
// This is a Rasmus Paulsen hack where there is direct access to the internal Quarternion and transformation
// .SECTION Caveats
// Whenever you add, subtract, or set points you must call Modified()
// on the vtkPoints object, or the transformation might not update.
// .SECTION see also
// vtkLinearTransform

#ifndef __vtkLandmarkTransformDirectQuart_h
#define __vtkLandmarkTransformDirectQuart_h

#include "vtkLinearTransform.h"

#define VTK_LANDMARK_RIGIDBODY 6
#define VTK_LANDMARK_SIMILARITY 7
#define VTK_LANDMARK_AFFINE 12

class vtkLandmarkTransformDirectQuart : public vtkLinearTransform
{
public:
  static vtkLandmarkTransformDirectQuart *New();

  vtkTypeMacro(vtkLandmarkTransformDirectQuart,vtkLinearTransform);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the source and target landmark sets. The two sets must have 
  // the same number of points.  If you add or change points in these objects,
  // you must call Modified() on them or the transformation might not update.
  void SetSourceLandmarks(vtkPoints *points);
  void SetTargetLandmarks(vtkPoints *points);
  vtkGetObjectMacro(SourceLandmarks, vtkPoints);
  vtkGetObjectMacro(TargetLandmarks, vtkPoints);

  vtkGetObjectMacro(Quarternion, double);
  vtkGetObjectMacro(Translation, double);


  // Description:
  // Set the number of degrees of freedom to constrain the solution to.
  // Rigidbody: rotation and Translation only.  
  // Similarity: rotation, Translation and isotropic scaling.
  // The default is similarity.
  vtkSetMacro(Mode,int);
  void SetModeToRigidBody() { this->SetMode(VTK_LANDMARK_RIGIDBODY); };
  void SetModeToSimilarity() { this->SetMode(VTK_LANDMARK_SIMILARITY); };
  void SetModeToAffine() { this->SetMode(VTK_LANDMARK_AFFINE); };
  vtkGetMacro(Mode,int);
  const char *GetModeAsString();

  // Description:
  // Invert the transformation.  This is done by switching the
  // source and target landmarks.
  void Inverse();

  // Description:
  // Get the MTime.
  vtkMTimeType GetMTime();

  // Description:
  // Make another transform of the same type.
  vtkAbstractTransform *MakeTransform();

protected:
  vtkLandmarkTransformDirectQuart();
  ~vtkLandmarkTransformDirectQuart();

  // Update the matrix from the quaternion.
  void InternalUpdate();

  // Description:
  // This method does no type checking, use DeepCopy instead.
  void InternalDeepCopy(vtkAbstractTransform *transform);

  vtkPoints* SourceLandmarks;
  vtkPoints* TargetLandmarks;

  int Mode;

  // Copy of internal Quarternion w x y z
  double Quarternion[4];

  // Copy of internal Translation (dx, dy, dz)
  double Translation[3];

private:
  vtkLandmarkTransformDirectQuart(const vtkLandmarkTransformDirectQuart&);  // Not implemented.
  void operator=(const vtkLandmarkTransformDirectQuart&);  // Not implemented.
};
 
//BTX
inline const char *vtkLandmarkTransformDirectQuart::GetModeAsString()
{
  switch (this->Mode)
    {
    case VTK_LANDMARK_RIGIDBODY:
      return "RigidBody";
    case VTK_LANDMARK_SIMILARITY:
      return "Similarity";
    case VTK_LANDMARK_AFFINE:
      return "Affine";
    default:
      return "Unrecognized";
    }
}
//ETX
#endif
