#ifndef __vtkRemoveNonManifoldCells_h
#define __vtkRemoveNonManifoldCells_h

#include <vtkPolyDataAlgorithm.h>

//! Remove triangles that are connected to a non-manifold edge
class vtkRemoveNonManifoldCells : public vtkPolyDataAlgorithm
{
public:
  static vtkRemoveNonManifoldCells *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkRemoveNonManifoldCells,vtkPolyDataAlgorithm);

  // Description:
  // Get the MTime of this object
  vtkMTimeType GetMTime();

  
protected:
  vtkRemoveNonManifoldCells();
 ~vtkRemoveNonManifoldCells();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkRemoveNonManifoldCells(const vtkRemoveNonManifoldCells&);  // Not implemented.
  void operator=(const vtkRemoveNonManifoldCells&);  // Not implemented.
};

#endif
