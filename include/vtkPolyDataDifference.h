#ifndef __vtkPolyDataDifference_h
#define __vtkPolyDataDifference_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkCellLocator.h"
#include "vtkpolydata.h"

class vtkPolyDataDifference : public vtkPolyDataAlgorithm
{
public:

  vtkTypeMacro(vtkPolyDataDifference,vtkPolyDataAlgorithm);
  static vtkPolyDataDifference *New();
  void PrintSelf(ostream& os, vtkIndent indent);


  void SetTargetData(vtkPolyData *pd); 
  vtkPolyData *GetTarget();

	vtkSetMacro(SignedDistance, bool);
	vtkGetMacro(SignedDistance, bool);

	vtkSetMacro(ExcludeEdgePoints, bool);
	vtkGetMacro(ExcludeEdgePoints, bool);

	vtkSetMacro(MaximumNumberOfLandmarks, int);
	vtkGetMacro(MaximumNumberOfLandmarks, int);

protected:
  vtkPolyDataDifference();
  ~vtkPolyDataDifference();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
   virtual int FillInputPortInformation(int port, vtkInformation *info);


  void CreateFeatureEdges(vtkPolyData *target);
 
  void DifferenceWithEdgePointsExcluded(vtkInformationVector **inputVector,  vtkInformationVector *outputVector);

  bool SignedDistance;

  int MaximumNumberOfLandmarks;

  // Should points falling on the edge of the target mesh be removed
  bool ExcludeEdgePoints;

  vtkPolyData *featureEdges;

  vtkCellLocator *featureLocator;

  // Tolerance for finding edges
  double EdgeFinderTolerance;

private:
  vtkPolyDataDifference(const vtkPolyDataDifference&) {};
  void operator=(const vtkPolyDataDifference&) {};

};

#endif
