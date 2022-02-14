#ifndef __vtkSignedDistanceTransformFilter_h
#define __vtkSignedDistanceTransformFilter_h

#include "vtkImageAlgorithm.h"
#include <vector>

class vtkPlane;
class vtkDoubleArray;

#define VTK_SDM_MODE_BRUTEFORCE 0
#define VTK_SDM_MODE_LOCALCELL 1
#define VTK_SDM_MODE_REALBRUTE 2


class vtkSignedDistanceTransformFilter : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkSignedDistanceTransformFilter,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSignedDistanceTransformFilter* New();

  // Description: 
  // Specify the spacing of the 3D sampling grid. If not set, a
  // reasonable guess will be made.
  vtkGetMacro(SampleSpacing,double);
  vtkSetMacro(SampleSpacing,double);

  	//! Insert a bounding plane, that are checked if points are inside or outside
	/** The normal of the plane should point into the region of interrest*/
	void InsertBoundingplane(vtkPlane *p);

  vtkSetClampMacro(SDMMode,int,
                   VTK_SDM_MODE_BRUTEFORCE,VTK_SDM_MODE_REALBRUTE);
  vtkGetMacro(SDMMode,int);

protected:
  vtkSignedDistanceTransformFilter();
  ~vtkSignedDistanceTransformFilter() {};

  virtual int RequestInformation (vtkInformation *,
                                  vtkInformationVector **,
                                  vtkInformationVector *);
  virtual int RequestData (vtkInformation *,
                           vtkInformationVector **,
                           vtkInformationVector *);

  virtual int FillInputPortInformation(int, vtkInformation*);


//  int NeighborhoodSize;
  double SampleSpacing;

  double bounds[6];

  double topleft[3];

  double bottomright[3];

  int dim[3];

private:
  vtkSignedDistanceTransformFilter(const vtkSignedDistanceTransformFilter&);  // Not implemented.
  void operator=(const vtkSignedDistanceTransformFilter&);  // Not implemented.

  vtkDoubleArray * AllocateVolume(vtkDataSet *input, vtkImageData *output, vtkInformation *outInfo, double Boundfactor);

	//! Calculate SDM by using a vtkCellLocator to find the closest point
	void ComputeSDMBruteforce(vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	//! Calculate SDM by finding distances in a band around all cells in the data set
	void ComputeSDMByCellProbing(vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	//! Calculate SDM by using a vtkCellLocator to find the closest point
	void ComputeSDMReallyBruteforce(vtkInformationVector** inputVector, vtkInformationVector* outputVector);

  	//! A set of planes that are checked to see if the points are outside
	/** The normals of the planes should point into the region of interrest*/
	std::vector<vtkPlane*> m_BoundingPlanes;

	//! The minimum bounding box size for a cell
	double m_MinBoundLength;

	//! Small class to store int points
	class XYZPoint
	{
	public:
		int x;
		int y;
		int z;

		XYZPoint(int xt, int yt, int zt) {x = xt; y = yt; z = zt;};
	};

	int SDMMode;

};

#endif

