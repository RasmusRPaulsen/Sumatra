// .NAME vtkPolyDataPCANormals - compute normals for polygonal mesh using a local PCA based method
// .SECTION Description
// vtkPolyDataPCANormals is a filter that computes point normals for a polygonal 
// mesh. A local PCA filter is used to determine the local normal.
// Only points with a well defined normal is returned
// Points will not be returned if: 
//  1. The points do not have enough neighbours to compute a PCA
//  2. The third Eigenvalue explains too much of the variation
//  3. The point is too far away from the estimated plane
//
// The returned scalar value per point is either the value of the 3rd Eigenvalue or the distance to the plane

#ifndef __vtkPolyDataPCANormals_h
#define __vtkPolyDataPCANormals_h

#include "vtkPolyDataAlgorithm.h"

#include <vector>

#define VTK_SCALAR_MODE_3RDEIGENVALUE 0
#define VTK_SCALAR_MODE_PLANEDIST 1

#define VTK_SEARCH_MODE_SEARCHRADIUS 0
#define VTK_SEARCH_MODE_NNEIGHBOURS 1

class  vtkPolyDataPCANormals : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPolyDataPCANormals,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	static vtkPolyDataPCANormals *New();

	// Description: 
	// The radius used when searching for neighbours
	vtkGetMacro(SearchRadius,double);
	vtkSetMacro(SearchRadius,double);

	// Description:
	// The number of points used to calculate a normal if the search is determined by number of neighbours
	vtkGetMacro(PointsPerNormals,int);
	vtkSetMacro(PointsPerNormals,int);

	// Description:
	// Threshold on the maximum distance from the point to the estimated plane
	vtkGetMacro(MaxPlaneDistance,double);
	vtkSetMacro(MaxPlaneDistance,double);


	// Description:
	// When the local PCA is done, the 3rd eigenvalue explains the amount of noise.
	// This sets a threshold on what is accepted
	vtkGetMacro(Max3EigenValue,double);
	vtkSetMacro(Max3EigenValue,double);

	// Description
	// Scalars are either the the 3rd eigenvalue or the current points distance to the PCA plane
	vtkSetClampMacro(ScalarMode,int,
		VTK_SCALAR_MODE_3RDEIGENVALUE,VTK_SCALAR_MODE_PLANEDIST);
	vtkGetMacro(ScalarMode,int);

	// Description
	// Search for neighbours is done using either using a search radius or a specific number of points to find
	vtkSetClampMacro(SearchMode,int,
		VTK_SEARCH_MODE_SEARCHRADIUS,VTK_SEARCH_MODE_NNEIGHBOURS);
	vtkGetMacro(SearchMode,int);


protected:
	vtkPolyDataPCANormals();
	~vtkPolyDataPCANormals() {};

	// Usual data generation method
	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	vtkPolyDataPCANormals(const vtkPolyDataPCANormals&);  // Not implemented.
	void operator=(const vtkPolyDataPCANormals&);  // Not implemented.

	//! Use a point locator to find the neighbouring points of all points
	/** Use the searchradius to find points */
	void CreateConnectivityUsingSearchRadius(vtkPolyData *input);

	//! Use a point locator to find the neighbouring points of all points
	void CreateConnectivityByPointsPerNormal(vtkPolyData *input);


	//! Estimate local normals using PCA
	void ComputeLocalPCA(vtkPolyData *input, vtkPolyData *output);
	
	//! Estimate local normals using PCA using a search radius to find neighbours
	void ComputeLocalPCAUsingSearchRadius(vtkPolyData *input, vtkPolyData *output);
	
	//! Estimate local normals using PCA using a fixed number of points per point
	void ComputeLocalPCAByPointsPerNormal(vtkPolyData *input, vtkPolyData *output);
	
	
	//! Locator search radius for finding connectivity
	double SearchRadius;

	//! The maximum percent the 3rd eigenvalue must explain of the local variation
	double Max3EigenValue;

	//! Keep track of neighbours of all points
	std::vector<std::vector<int> > m_Connectivity;

	int ScalarMode;

	// The number of points used to calculate a normal if the search is determined by number of neighbours
	int PointsPerNormals;

	int SearchMode;

	//! Maximum allowed distance to estimated plane
	double MaxPlaneDistance;
};

#endif
