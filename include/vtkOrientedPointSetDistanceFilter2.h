#ifndef __vtkOrientedPointSetDistanceFilter2_h
#define __vtkOrientedPointSetDistanceFilter2_h

#include "vtkImageAlgorithm.h"

class vtkDoubleArray;
class vtkImageData;
class vtkIntArray;
class vtkPointLocator;
class vtkCellLocator;

#define VTK_ORIENTEDDISTANCE_SIMPLE 0
#define VTK_ORIENTEDDISTANCE_PROJECTED 1
#define VTK_ORIENTEDDISTANCE_PROJECTED_MEDIAN 2
#define VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE 3
#define VTK_ORIENTEDDISTANCE_AVERAGE 4
#define VTK_ORIENTEDDISTANCE_MEDIAN 5
#define VTK_UNSIGNEDDISTANCE_MEDIAN 6
#define VTK_UNSIGNEDDISTANCE_SURFACE 7


#define VTK_ORIENTEDDISTANCE_FULLVOLUME 0
#define VTK_ORIENTEDDISTANCE_BAND 1
#define VTK_ORIENTEDDISTANCE_COMBINED_FULLBAND 2
#define VTK_UNSIGNEDDISTANCE_FULLVOLUME 3

#define VTK_SEARCH_MODE_SEARCHRADIUS 0
#define VTK_SEARCH_MODE_NNEIGHBOURS 1


//! Create a volume and for each voxel calculate the signed distance to the input polydata
/** The input is considered a point set, where each point also has an associated normals.
Two different distance measures can be used. Either the simple Euclidian distance between the voxel
and the closest point or a "plane distance", where the normal of the point is taken into account.
Compared to vtkOrientedPointSetDistanceFilter this class returns a weight volume
and can be used to update an existing distance field */
class vtkOrientedPointSetDistanceFilter2 : public vtkImageAlgorithm
{
public:
	vtkTypeMacro(vtkOrientedPointSetDistanceFilter2,vtkImageAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	// Construct with NeighborhoodSize=20.
	static vtkOrientedPointSetDistanceFilter2* New();

	// Description: 
	// The radius used when searching for neighbours
	vtkGetMacro(SearchRadius,double);
	vtkSetMacro(SearchRadius,double);

	// Description: 
	//! The value used in normalising the distance weight
	/** It is also determines how far out should be marked in the band distance update */
	vtkGetMacro(WeightValueHigh,double);
	vtkSetMacro(WeightValueHigh,double);

	// Description
	// Search for neighbours is done using either using a search radius or a specific number of points to find
	vtkSetClampMacro(SearchMode,int,
		VTK_SEARCH_MODE_SEARCHRADIUS,VTK_SEARCH_MODE_NNEIGHBOURS);
	vtkGetMacro(SearchMode,int);


	// Description: 
	// How many distances should be used in the distance calculation
	vtkGetMacro(NumberOfDistances,int);
	vtkSetMacro(NumberOfDistances,int);

	// Description: 
	// How many threads can be used. -1 single threaded. 0 automatically determine number of threads (#cores-1)
	vtkGetMacro(NumThreads,int);
	vtkSetMacro(NumThreads,int);

	//! Distance mode using when creating the distance map
	/** 0 = Nearest point Euclidean distance (currently equal to 4), 
		1 = Nearest point projected distance (currently equal to 3), 
		2 = Median of nearest points projected distance (L1 norm),
		3 = Average of nearest points projected distance (L2 norm)
		4 = Average of distances to nearest points
		5 = Median of distances to nearest points
		6 = Median of unsigned distance to nearest points
		7 = Unsigned distance to closest surface */
	vtkSetClampMacro(DistanceMode,int,
		VTK_ORIENTEDDISTANCE_SIMPLE, VTK_UNSIGNEDDISTANCE_SURFACE);
	vtkGetMacro(DistanceMode,int);

	//! The modes of operations
	/** Either the full volume is computed or an existing volume is updated in a band around the points.
		In the third mode a new volume is created but only a band is calculated*/
	vtkSetClampMacro(ComputeMode,int,
		VTK_ORIENTEDDISTANCE_FULLVOLUME,VTK_UNSIGNEDDISTANCE_FULLVOLUME);
	vtkGetMacro(ComputeMode,int);

	// Description:
	// Should a reference volume be created
	// Each voxel in the reference volume contains the id of the point that
	// are used to create the distance in the distance volume
	// This can be used to create weighting volumes for further processing
	vtkSetMacro(CreateWeightVolume,int);
	vtkGetMacro(CreateWeightVolume,int);
	vtkBooleanMacro(CreateWeightVolume,int);

	//! Get pointer to reference volume
	vtkGetMacro(WeightVolume,vtkImageData*);

	//! Set pointer to existing distance map
	vtkSetMacro(InputDistances,vtkImageData*);

	//! How many voxels should the min length side have when constructing volume
	vtkGetMacro(MinSideVoxels,int);
	vtkSetMacro(MinSideVoxels,int);

	// Description: 
	//! How many voxels should be on each side of the point cloud in the volume
	vtkGetMacro(PadVoxels,int);
	vtkSetMacro(PadVoxels,int);


protected:
	vtkOrientedPointSetDistanceFilter2();
	~vtkOrientedPointSetDistanceFilter2();

  virtual int RequestInformation (vtkInformation *,
                                  vtkInformationVector **,
                                  vtkInformationVector *);
  virtual int RequestData (vtkInformation *,
                           vtkInformationVector **,
                           vtkInformationVector *);

  virtual int FillInputPortInformation(int, vtkInformation*);

	//! Calculate the distance from a point to a set of closest neighbours
	/** Different methods can be used */
	bool CalculateDistanceAtPosition( vtkPointLocator * locator, double * point, vtkDataSet * input, vtkDataArray * normals, vtkDoubleArray * newScalars, int offset, int LocalSearchMode );

	//! Calculate the distance from a point to the surface
	bool CalculateDistanceToSurfaceAtPosition(vtkCellLocator* locator, double* point, vtkDataSet* input, vtkDataArray* normals, vtkDoubleArray* newScalars, int offset);

	//! Compute distances for the full distance field
	void ComputeDistancesInFullVolume(vtkDataSet *input, int dim[3],
		vtkDoubleArray *newScalars, double topleft[3],
		vtkDataArray *normals);
	
	//! Compute distances for the full distance field. Multi-threaded
	void ComputeDistancesInFullVolumeMultiThreaded(vtkDataSet *input, int dim[3], vtkDoubleArray *newScalars, double topleft[3], vtkDataArray *normals);

	
	//! Compute distances in the distance field in a band around input points. Interpolate values other places using upsampled input
	void ComputeDistancesInBandWithInput(vtkDataSet *input,  int dim[3], vtkImageData* output, double topleft[3], vtkDataArray *normals);

	//! Compute distances in the distance field in a band around input points. Interpolate values other places using upsampled input. Multithreaded
	void ComputeDistancesInBandWithInputMultiThreaded(vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDataArray *normals);

	//! Compute distances in the distance field in a band around input points. Values outside bands are set to 0
	void ComputeDistancesInBand(vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDataArray *normals);
	
	//! Mark voxels to be updated
	void BandMarkVolumeInBand(vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDoubleArray* MarkVolumeData);

	//!  int NeighborhoodSize;
	double SampleSpacing;

private:
	vtkOrientedPointSetDistanceFilter2(const vtkOrientedPointSetDistanceFilter2&);  // Not implemented.
	void operator=(const vtkOrientedPointSetDistanceFilter2&);  // Not implemented.

	//! Create a new output volume - a little smarter
	void CreateNewVolumeFromPolyData(vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	void CreateNewVolumeFromInputVolume(vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	// Each voxel in the weight volume contains a weight
	vtkImageData	*WeightVolume;

	vtkDoubleArray		*WeightVolumeData;

	vtkImageData *InputDistances;

	// This factor determines the scaling of the size of the voxels
	double SampleFactor;

	int DistanceMode;

	//! Should reference volume be computed
	int CreateWeightVolume;

	//! How many distances should be used in the median computation
	int NumberOfDistances;

	//! The value used in normalising the distance weight
	double WeightValueHigh;

	int ComputeMode;

	//! How many percent compared to the polydata bounding box should the volume be enlarged
	double PadVoxels;

	//! How many voxels should the min length side have when constructing volume
	int MinSideVoxels;

	int SearchMode;

	double SearchRadius;

	// How many threads can be used. -1 single threaded. 0 automatically determine number of threads (#cores-1)
	int NumThreads;
};

#endif

