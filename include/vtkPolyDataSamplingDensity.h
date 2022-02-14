// .NAME vtkPolyDataSamplingDensity - Compute the local point density and return it as the scalar values
// .SECTION Description
// vtkPolyDataSamplingDensity  

#ifndef __vtkPolyDataSamplingDensity_h
#define __vtkPolyDataSamplingDensity_h

#include "vtkPolyDataAlgorithm.h"

#include <vector>

class  vtkPolyDataSamplingDensity : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPolyDataSamplingDensity,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	static vtkPolyDataSamplingDensity *New();

	// Description: 
	// The radius used when searching for neighbours
	vtkGetMacro(SearchRadius,double);
	vtkSetMacro(SearchRadius,double);


protected:
	vtkPolyDataSamplingDensity();
	~vtkPolyDataSamplingDensity() {};

	// Usual data generation method
	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	vtkPolyDataSamplingDensity(const vtkPolyDataSamplingDensity&);  // Not implemented.
	void operator=(const vtkPolyDataSamplingDensity&);  // Not implemented.

	//! Estimate local density
	void ComputeLocalDensity(vtkPolyData *input, vtkPolyData *output);

	//! Locator search radius for finding connectivity
	double SearchRadius;
};

#endif
