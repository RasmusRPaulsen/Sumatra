#ifndef _vtkExtMisc_h_
#define _vtkExtMisc_h_

#include <iostream>
#include <string>

class vtkPolyData;
class vtkMatrix4x4;
class vtkPlaneSource;
class vtkPoints;
class vtkIdList;
class vtkCell;
class vtkPolyDataReader;
class vtkSTLReader;
class vtkOBJReader;
class vtkPlane;
class vtkRenderWindow;
class vtkDoubleArray;
class vtkFloatArray;

//! General functions that can be used together with VTK
class vtkExtMisc
{
public:

	//! Generate a transformation matrix that will transform a shape to be aligned with the line
	/** \param theta the rotation around the axis */
	static void TransformMatrixFromLine(double *p1, double *p2, double theta, vtkMatrix4x4* mat);

	//! Print info about polydata
	static void PrintPolyInfo(vtkPolyData* pd);

	//! Write Polydata to a VTK file
	/** Returns false in case of errors */
	static bool WritePDVTK(vtkPolyData *pd, const std::string& fname);
	
	//! Write points as spheres with a given radius
	/** if radius == 0 an estimate is made of the best radius */
	static bool WritePDFPointsAsSpheres(vtkPoints *points, double radius, const std::string& fname);

	//! Write Polydata to a STL file
	/** Returns false in case of errors */
	static bool WritePDSTL(vtkPolyData *pd, const std::string& fname, bool binary = false);

	//! Write Polydata to a STL file
	/** Returns false in case of errors */
	static bool WritePDSTLWithNormals(vtkPolyData *pd, const std::string& fname, bool binary = false);

	//! Write Polydata to a OBJ file
	/** Notice that this method uses the render to export the whole scene
	    Returns false in case of errors */
	static bool WritePDOBJ(vtkRenderWindow *ren, const std::string& fname);

	//! Generate normals and write Polydata to a VTK file
	static void WritePDWithNormals(vtkPolyData *pd, const std::string& fname);

	//! Write a vtkPoints as a polydata with vertex structures
	static void WritePointsAsVertices(vtkPoints *pts, const std::string& fname);

	//! Write a point cloud as a fitted blob
	static void WritePointsAsPCAFittedBlob(vtkPoints *pts, const std::string &fname);

	//! Allocate matrix
	static double** NewMatrix(int rows, int cols);

	//! Delete matrix
	static void DeleteMatrix(double **m);

	//! Set matrix to zero
	static void ZeroMatrix(double **m, int rows, int cols);

	//! Multiply matrix
	static void MatrixMultiply(double **a, double **b, double **c,
				  int arows, int acols, 
				  int brows, int bcols);

	//! Create large covariance matrix
	static void LargeCovarianceMatrix(double **a, double **c,
				  int arows, int acols);

	//! Create small covariance matrix
	static void SmallCovarianceMatrix(double **a, double **c,
				  int arows, int acols);

	//! Calculates the "small" covariance for the points pts
	/** If premean, precov, prevNumberOfPoints are not NULL
	    then they are assumed to contain the mean, covariance and number of points
		for the first part of the data set
		Then newmean, newcov and NumberOfPoints will be the mean, covariance and number of points
		for the entire data set					*/
	static void IterativeCovarianceMatrix(vtkPoints *pts, 
				  double *prevMean, double **prevCov, int prevNumberOfPoints, 
				  double *newMean, double ** newCov, int *numberOfPoints);

	//! Transpose matrix
	static void MatrixTranspose(double **a, double **b, int rows, int cols);

	//! Write matrix as ASCII
	static void WriteMatrixASCII(double **m, int rows, int cols, const std::string& fname);
	
	//! Write matrix as a vtk format file
	static void WriteMatrixVTKFormat(double **m, int rows, int cols, const std::string& fname, bool binary);

	//! Read matrix from a VTK format file
	static double **ReadMatrixVTKFormat(int &rows, int &cols, const std::string& fname);

	//! Write an array to a simple text file
	static void WriteFloatArray(vtkFloatArray *a, const std::string &fname);

	//! Read a simple text array
	static bool ReadDoubleArray(const std::string& fname, vtkDoubleArray *array);

	//! Allocate vector
	static double* NewVector(int length);

	//! Delete vector
	static void DeleteVector(double *v);

	//! Fit a plane to a collection of points.
	/** returned as a planeSource to be able to use the spanning points
	    \param size used to scale the plane. If size = 1 the plane size is equal to sqrt(eigenvalues).*/
	static void FitPlaneToPointsLSQ(vtkPoints *pts, vtkPlaneSource *plane, double size);

	//! Compute the spanning vectors of a point cloud
	/** Eigenvalues and the mean of the point cloud are also returned.
	    NOTE: Currently the eigenvectors are returnes "transposed", so 
	    to get the elements of the first vector use (v0[0], v1[0], v2[0]).*/
	static void ComputeEigenVectorsOfPointCloud(vtkPoints *pts, double* mean, double *v0, double *v1, double *v2,
												 double* evals);

	//! Compute least-squares intersection points for a set of lines
	/** Given a set of lines defined by start and end points the point with the 
	    least sum of squared distances to the lines is returned */
	static void ComputeLeastSquaresIntersectionPoint(vtkPoints *startPoints, vtkPoints *endPoints, double *p);
	
	//! Compute least-squares intersection points for a set of lines
	/** Use RANSAC to remove outliers */
	static void ComputeLeastSquaresIntersectionPointWithRANSAC(vtkPoints *startPoints, vtkPoints *endPoints, double *p);
	
	//! Save a plane as corners												 
	static void SavePlaneAsCorners(vtkPlaneSource *plane, const std::string& fname);

	//! Calculate a planes sidelength and the ABCD representation
	/** based on three corner points.Origin (po), span point 1 (p1) and span point 2 (p2)
	    ABCD are normalised so the vector ABCD have norm 1*/
	static void Plane2Parameters(double *po, double *p1, double *p2, double &SideLength, double &A, double &B, double &C, double &D);

	//! Calculate the list of connected neighbours to a given vertex
	/** all connected neighbours with a squared Euclidian distance less than Neighborhoodise2 is included
	    requires that mesh->BuildLinks() has been called*/
	static void FindVertexNeighbours(vtkPolyData* mesh, int VertexID, vtkIdList *NeighbourIDs, double NeighborhoodSize2);

	//! Read a polydata and check if the file exists
	static vtkPolyDataReader* SafeReadPolyData(const std::string& fname);

	//! Read a polydata as a STL file and check if the file exists
	static vtkSTLReader* SafeReadPolyDataSTL(const std::string& fname);

	//! Read a polydata as a OBJ file and check if the file exists
	static vtkOBJReader* SafeReadPolyDataOBJ(const std::string& fname);

	//! Calculate the center of mass of a mesh
	static void CenterOfMass(vtkPolyData* pd, double *CM);

	//! Calculate the center of mass of a mesh
	static void CenterOfMass(vtkPoints* pd, double *CM);

	//! Calculate the center of mass of a mesh, but only for the vertices in the list
	static void CenterOfMass(vtkPolyData* pd, vtkIdList* pids, double *CM);

	//! Translate the center of mass of a mesh to origo
	static void TranslateCenterOfMassToOrigo(vtkPolyData *inPD, vtkPolyData *outPD);

	//! Translate the center of mass to origo and scale the mesh so the largest side length is equal to size
	/** The center of mass is returned together with the scale factor used to scale the mesh */
	static void NormaliseSize(vtkPolyData *inPD, vtkPolyData *outPD, double size, double *CM, double &scale);

	//! Calculate the area of a planar polygon
	static double PolygonArea(vtkPoints *poly);

	//! Calculate the area of a planar polygon using a 2D projection technique
	static double PolygonAreaByProjection(vtkPoints *poly);

	static double PolygonAreaByProjection(vtkCell *poly);

	static double PolygonAreaByProjection(vtkPoints *poly, double *Normal);

	//! For all points calculate the distance to the center of mass
	/** return the average and sdev.
	   Can be used to see how close a cut approximates a complete circle*/
	static void DistanceFromCMStats(vtkPoints *poly, double &DistMean, double &DistSDev);

	//! Make a 2D projection and write it to a file
	static void ProjectPolygonTo2DWriteToFile(vtkPoints *poly, const std::string fname);

	//! Find the intersection (if any) between a plane and a polyline
	/** the plane is defined as an origin and a normal*/
	static bool InterSectionBetweenPolyLineAndPlane(vtkPolyData* Pline, double *pNorm, double *pOrg,
																int& lowPointId, double &IntersectT, double *iPoint);


	//! Move distance dist along a polyline. The start point is defined using lowPointId and IntersecT
	/** return false if the distance extends beyond the polyline*/
	static bool MoveDistanceAlongPolyLine(vtkPolyData* Pline, int lowPointId, double t, double dist,
		int& NewLowPointId, double &NewT, double *NewP);


	//! Given a point, porg, and a normal direction, pnoraml, in that point, determine if another point, p, is on the inside or the outside
	/** \return positive number if point is placed in the normal direction else negative*/
	static double PointInsideOrOutside(double *porg, double *pnormal, double *p);

	//! Cut a poly data by a plane and calculate the area of the resulting cut
	/** Assumes that the resulting cut is equivalent to a closed circle */
	static double CutAndCalculateArea(vtkPolyData *pd, vtkPlane *plane);

	//! Calculate the area of the polydata
	static double SurfaceArea(vtkPolyData *pd);

	//! Find the closest point (cp) on the surface from point x and return the distance
	/** Assumes a triangulated surface */
	static double FindClosestPointOnSurface(vtkPolyData *pd, double *x, double *cp);

	//! Returns some computed surface values in string format
	/** \param sep separator char to use
	    \param spa space char to use.*/
	static std::string GetSurfaceValues(vtkPolyData *pd, const std::string& sep, const std::string& spa = "");

	//! Check if the surface is consistent (non-manifold edges etc)
	/** Returns true if ok. msg contains error message */
	static bool CheckSurfaceConsistency(vtkPolyData *pd, std::string& msg);

	//! Can read surface data in several different formats
	static vtkPolyData *MultiReadSurface(const std::string &fname);

	//! Can write surface data in several different formats
	static bool MultiWriteSurface(const std::string& fname, vtkPolyData* pd, bool WriteNormals = false, bool WriteScalars = false);
};

#endif
