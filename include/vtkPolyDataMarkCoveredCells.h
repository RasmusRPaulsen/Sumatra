#ifndef __vtkPolyDataMarkCoveredCells_h
#define __vtkPolyDataMarkCoveredCells_h

#include <vtkPolyDataAlgorithm.h>

class vtkPolyData;

//! Mark points that are near points in the source polydata
/** Can be used to later remove points that are far away from a source mesh */
class vtkPolyDataMarkCoveredCells : public vtkPolyDataAlgorithm
{
public:
  static vtkPolyDataMarkCoveredCells *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkPolyDataMarkCoveredCells,vtkPolyDataAlgorithm);

  // Description:
  // Get the MTime of this object
  vtkMTimeType GetMTime();
  
  // Description: 
  vtkSetMacro(Source,vtkPolyData*);
  vtkGetMacro(Source,vtkPolyData*);

  // Description: 
  // The radius used to mark points
  vtkSetMacro(MarkRadius,double);
  vtkGetMacro(MarkRadius,double);

protected:
  vtkPolyDataMarkCoveredCells();
 ~vtkPolyDataMarkCoveredCells();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkPolyDataMarkCoveredCells(const vtkPolyDataMarkCoveredCells&);  // Not implemented.
  void operator=(const vtkPolyDataMarkCoveredCells&);  // Not implemented.

  double MarkRadius;

  vtkPolyData *Source;
};

#endif
