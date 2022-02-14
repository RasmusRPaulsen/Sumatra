// .NAME vtkPolyDataResamplePoints - samples a given number of points from a point cloud
// .SECTION Description
// vtkPolyDataResamplePoints is a filter that can reduce the number of points in a point cloud
// it simply copies the required number of points by sampling a fixed interval

#ifndef __vtkPolyDataResamplePoints_h
#define __vtkPolyDataResamplePoints_h

#include "vtkPolyDataAlgorithm.h"

#include <vector>

class  vtkPolyDataResamplePoints : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPolyDataResamplePoints,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	static vtkPolyDataResamplePoints *New();

	// Description:
	// The number of point in the output cloud
	vtkGetMacro(NumberOfOutputPoints,int);
	vtkSetMacro(NumberOfOutputPoints,int);


protected:
	vtkPolyDataResamplePoints();
	~vtkPolyDataResamplePoints() {};

	// Usual data generation method
	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	vtkPolyDataResamplePoints(const vtkPolyDataResamplePoints&);  // Not implemented.
	void operator=(const vtkPolyDataResamplePoints&);  // Not implemented.

	//! Number of output points
	int NumberOfOutputPoints;
};

#endif
