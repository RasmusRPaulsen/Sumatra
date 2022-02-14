#ifndef __vtkRemoveUnusedPolyDataPoints_h
#define __vtkRemoveUnusedPolyDataPoints_h

#include "vtkPolyDataAlgorithm.h"

class vtkRemoveUnusedPolyDataPoints : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkRemoveUnusedPolyDataPoints,vtkPolyDataAlgorithm);
  static vtkRemoveUnusedPolyDataPoints *New();
  void PrintSelf(ostream& os, vtkIndent indent);


protected:
	vtkRemoveUnusedPolyDataPoints();
  ~vtkRemoveUnusedPolyDataPoints();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkRemoveUnusedPolyDataPoints(const vtkRemoveUnusedPolyDataPoints&);
  void operator=(const vtkRemoveUnusedPolyDataPoints&);
};

#endif


