#ifndef __vtkCropPolyData_h
#define __vtkCropPolyData_h

#include <vtkPolyDataAlgorithm.h>

#include <vtkPolyData.h>
//#include <vector>

//! Graph based cropping.
/** Removes points that have a scalar value of 0 and that are connected to the bounding box. */
class vtkCropPolyData : public vtkPolyDataAlgorithm
{
public:
  static vtkCropPolyData *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkCropPolyData,vtkPolyDataAlgorithm);

	// Set boundary that is used for cropping
	void SetBoundary(double boundary[6]);

	vtkSetMacro(Tolerance,double);

  // Description:
  // Get the MTime of this object
	vtkMTimeType GetMTime();
 

protected:
  vtkCropPolyData();
 ~vtkCropPolyData();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkCropPolyData(const vtkCropPolyData&);  // Not implemented.
  void operator=(const vtkCropPolyData&);  // Not implemented.

  // Find points on the bounding box
  void LocateEdgePoints(vtkPolyData *input, vtkPolyData *output);

  void MarkConnectedPoints(vtkPolyData *output);

  // Status of each point. 0 - keep it, 1 - throw it away
//  std::vector<int> PointStatus;

  double Boundary[6];

  double Tolerance;
};

#endif
