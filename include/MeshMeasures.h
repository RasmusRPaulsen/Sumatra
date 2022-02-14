#ifndef _meshmeasures_h_
#define _meshmeasures_h_

#include <vector>
#include <vtkType.h>

class vtkPolyData;
class vtkCell;
class vtkDataSet;

//! Functions to measure properties of a mesh
class CMeshMeasures  
{
public:
	//! Calculate the average and sdev of edge lengths in a mesh
	static void MeshRegularityBasedOnEdgeLengths(vtkPolyData *mesh, double &DistMean, double &DistSdev);

	//! Calculate the average of the minimum angles in the mesh
	/** The minimum angle is the smallest angle in a triangle.*/
	static void MeshRegularityBasedOnMinimumTriangleAngle(vtkPolyData *mesh, double &MinAngleMean);

	//! Calculate all the minimum angles in the mesh
	/** The minimum angle is the smallest angle in a triangle.*/
	static void AllMinimumAngles(vtkPolyData *mesh, std::vector<double> &MinAngles);

	//! Examine if all points of cell1 are also in cell2
	/** returns true if if all points of cell1 are also in cell2.
		Also check for mirrored cells.*/
	static bool IsCellContainedInOtherCell(vtkCell *cell1, vtkCell *cell2);

	//! Examine if all points of cell1 are also in cell2
	/** returns true if if all points of cell1 are also in cell2.
		Also check for mirrored cells.*/
	static bool IsCellContainedInOtherCell(vtkIdType npts1, const vtkIdType *pts1, vtkIdType npts2, const vtkIdType *pts2);

	//! Calculate the number of cells contained in other cells in the mesh
	static int TotalNumberOfCellsContainedInOtherCells(vtkPolyData *pd);

	//! Compute statistics on the closest neighbour point in point sests
	static bool NearestNeighbourStatistics(vtkDataSet *pd, double &minDist, double &maxDist, double &meanDist, double &sdvDist, double &medianDist, double &frac05, double &frac95);

	//! Calculate the three sidelengths of a triangle and return false if one is 0
	static bool TriangleSideLengths(vtkPolyData *mesh, int cellId, double &l1, double &l2, double &l3);

	//! Calculate the three angles of an triangle and alse the minimum angle
	static bool TriangleAngles(vtkPolyData *mesh, int cellId, double &A, double &B, double &C, double &minA);
};

#endif
