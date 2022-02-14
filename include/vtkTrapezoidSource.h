
#ifndef __vtkTrapezoidSource_h
#define __vtkTrapezoidSource_h

#include "vtkPolyDataAlgorithm.h"

//! vtkTrapezoidSource creates a trapezoid centered at origin. 
/** The trapezoi is represented with four-sided polygons. */
class vtkTrapezoidSource : public vtkPolyDataAlgorithm 
{
public:
  static vtkTrapezoidSource *New();
  vtkTypeMacro(vtkTrapezoidSource,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the length of the trapezoid in the x-direction.
  vtkSetClampMacro(XLength,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(XLength,double);

  // Description:
  // Set the length of the trapezoid in the y-direction.
  vtkSetClampMacro(YLength,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(YLength,double);

  // Description:
  // Set the length of the trapezoid in the z-direction.
  vtkSetClampMacro(ZLength,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(ZLength,double);

  // Description:
  // Set the Alpha angle in the trapezoid
  vtkSetClampMacro(Alpha,double,-360,360);
  vtkGetMacro(Alpha,double);
  
   // Description:
  // Set the Beta angle in the trapezoid
  vtkSetClampMacro(Beta,double,-360,360);
  vtkGetMacro(Beta,double);

  // Description:
  // Set the center of the trapezoid.
  vtkSetVector3Macro(Center,double);
  vtkGetVectorMacro(Center,double,3);

  // Description:
  // Turn on/off whether to cap trapezoid with polygons.
  vtkSetMacro(Capping,int);
  vtkGetMacro(Capping,int);
  vtkBooleanMacro(Capping,int);

protected:
  vtkTrapezoidSource(double xL=1.0, double yL=1.0, double zL=1.0, double A = 90, double B = 90, int Cap = 1);
  ~vtkTrapezoidSource() {};

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  double XLength;
  double YLength;
  double ZLength;
  double Alpha;
  double Beta;
  double Center[3];
  int Capping;

private:
  vtkTrapezoidSource(const vtkTrapezoidSource&);  // Not implemented.
  void operator=(const vtkTrapezoidSource&);  // Not implemented.
};

#endif
