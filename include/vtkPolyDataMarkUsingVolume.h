#ifndef __vtkPolyDataMarkUsingVolume_h
#define __vtkPolyDataMarkUsingVolume_h

#include <vtkPolyDataAlgorithm.h>

class vtkPolyData;
class vtkImageData;

//! Mark a polydata based on the values of a voxel volume
/** Each point in the polydata is assigned the interpolated voxel value */
class vtkPolyDataMarkUsingVolume : public vtkPolyDataAlgorithm
{
public:
  static vtkPolyDataMarkUsingVolume *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkPolyDataMarkUsingVolume,vtkPolyDataAlgorithm);

  // Description:
  // Get the MTime of this object
  vtkMTimeType GetMTime();
 
  // Set the volume to use
  vtkSetMacro(Volume,vtkImageData*);

protected:
  vtkPolyDataMarkUsingVolume();
 ~vtkPolyDataMarkUsingVolume();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkPolyDataMarkUsingVolume(const vtkPolyDataMarkUsingVolume&);  // Not implemented.
  void operator=(const vtkPolyDataMarkUsingVolume&);  // Not implemented.

  vtkImageData *Volume;
};

#endif
