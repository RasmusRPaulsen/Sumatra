#ifndef __vtkPolyDataRandomResamplePoint_h
#define __vtkPolyDataRandomResamplePoint_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkCellLocator.h"

//! Resample a surface using random points
/** All connectivity is killed */
class vtkPolyDataRandomResamplePoint : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPolyDataRandomResamplePoint,vtkPolyDataAlgorithm);
  static vtkPolyDataRandomResamplePoint *New();
  void PrintSelf(ostream& os, vtkIndent indent);
  

  	vtkSetMacro(ResamplePoints, int);
	vtkGetMacro(ResamplePoints, int);

	// Description:
	vtkSetMacro(CreatePointNormals,int);
	vtkGetMacro(CreatePointNormals,int);
	vtkBooleanMacro(CreatePointNormals,int);

protected:
  vtkPolyDataRandomResamplePoint();
  ~vtkPolyDataRandomResamplePoint();

   virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  //! Number of points to sample
	int ResamplePoints;

	//! border to add around bounding box when sampling
	double AddBorder;

	int CreatePointNormals;

private:
  vtkPolyDataRandomResamplePoint(const vtkPolyDataRandomResamplePoint&) {};
  void operator=(const vtkPolyDataRandomResamplePoint&) {};
};


#endif
