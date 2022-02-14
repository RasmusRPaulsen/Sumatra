#ifndef __vtkOrientedPointSetDistanceFilter_h
#define __vtkOrientedPointSetDistanceFilter_h

#include "vtkImageAlgorithm.h"

class vtkDoubleArray;
class vtkImageData;
class vtkIntArray;


#define VTK_ORIENTEDDISTANCE_SIMPLE 0
#define VTK_ORIENTEDDISTANCE_PROJECTED 1
#define VTK_ORIENTEDDISTANCE_PROJECTED_MEDIAN 2
#define VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE 3

//! Create a volume and for each voxel calculate the signed distance to the input polydata
/** The input is considered a point set, where each point also has an associated normals.
Two different distance measures can be used. Either the simple Euclidian distance between the voxel
and the closest point or a "plane distance", where the normal of the point is taken into account.*/
class vtkOrientedPointSetDistanceFilter : public vtkImageAlgorithm
{
public:
	vtkTypeMacro(vtkOrientedPointSetDistanceFilter,vtkImageAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	// Construct with NeighborhoodSize=20.
	static vtkOrientedPointSetDistanceFilter* New();

	// Description: 
	// Specify the spacing of the 3D sampling grid. If not set, a
	// reasonable guess will be made.
	vtkGetMacro(SampleSpacing,double);
	vtkSetMacro(SampleSpacing,double);

	// Description: 
	// Specify the spacing of the 3D sampling grid. 
	// It scales the guess. The higher value, the finer grid
	vtkGetMacro(SampleFactor,double);
	vtkSetMacro(SampleFactor,double);

	// Description: 
	// How many distances should be used in the distance calculation
	vtkGetMacro(NumberOfDistances,int);
	vtkSetMacro(NumberOfDistances,int);

	//! Distance mode using when creating the distance map
	/** 0 = Nearest point Euclidean distance, 
		1 = Nearest point projected distance, 
		2 = Median of nearest points projected distance (L1 norm),
		3 = Average of nearest points projected distance (L2 norm)*/
	vtkSetClampMacro(DistanceMode,int,
		VTK_ORIENTEDDISTANCE_SIMPLE,VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE);
	vtkGetMacro(DistanceMode,int);

	// Description:
	// Should a reference volume be created
	// Each voxel in the reference volume contains the id of the point that
	// are used to create the distance in the distance volume
	// This can be used to create weighting volumes for further processing
	vtkSetMacro(CreateReferenceVolume,int);
	vtkGetMacro(CreateReferenceVolume,int);
	vtkBooleanMacro(CreateReferenceVolume,int);

	//! Get pointer to reference volume
	vtkGetMacro(ReferenceVolume,vtkImageData*);

protected:
	vtkOrientedPointSetDistanceFilter();
	~vtkOrientedPointSetDistanceFilter();

  virtual int RequestInformation (vtkInformation *,
                                  vtkInformationVector **,
                                  vtkInformationVector *);
  virtual int RequestData (vtkInformation *,
                           vtkInformationVector **,
                           vtkInformationVector *);

  virtual int FillInputPortInformation(int, vtkInformation*);


	//! Use the nearest neighbour found by the locator
	void UseSimpleCloseNeighbour(vtkDataSet *input, int dim[3],
		vtkDoubleArray *newScalars, double topleft[3],
		vtkDataArray *normals);

	//! Use the nearest neighbour found using a normal-projection method
	void UseProjectedCloseNeighbour(vtkDataSet *input, int dim[3],
		vtkDoubleArray *newScalars, double topleft[3],
		vtkDataArray *normals);

	//! Use the nearest neighbour found using a normal-projection method
	/** Furthermore, the median of distances is taken */
	void UseProjectedCloseNeighbourMedian(vtkDataSet *input, int dim[3],
		vtkDoubleArray *newScalars, double topleft[3],
		vtkDataArray *normals);

	//! Use the nearest neighbour found using a normal-projection method
	/** Furthermore, the average of distances is taken */
	void UseProjectedCloseNeighbourAverage(vtkDataSet *input, int dim[3],
		vtkDoubleArray *newScalars, double topleft[3],
		vtkDataArray *normals);

	//  int NeighborhoodSize;
	double SampleSpacing;
private:
	vtkOrientedPointSetDistanceFilter(const vtkOrientedPointSetDistanceFilter&);  // Not implemented.
	void operator=(const vtkOrientedPointSetDistanceFilter&);  // Not implemented.

	// Each voxel in the reference volume contains the id of the point that
	// are used to create the distance in the distance volume
	vtkImageData	*ReferenceVolume;

	vtkIntArray		*ReferenceVolumeData;

	// This factor determines the scaling of the size of the voxels
	double SampleFactor;

	int DistanceMode;

	//! Should reference volume be computed
	int CreateReferenceVolume;

	//! How many distances should be used in the median computation
	int NumberOfDistances;
};

#endif

