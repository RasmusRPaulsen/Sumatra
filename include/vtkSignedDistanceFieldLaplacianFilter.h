#ifndef __vtkSignedDistanceFieldLaplacianFilter_h
#define __vtkSignedDistanceFieldLaplacianFilter_h

#include "vtkSimpleImageToImageFilter.h"
#include <vector>

class vtkDoubleArray;

// Performs smoothing of the gradient of a signed distance field
// The approach is simalar to Markov Random Field optimisation
class vtkSignedDistanceFieldLaplacianFilter : public vtkSimpleImageToImageFilter
{
public:
	static vtkSignedDistanceFieldLaplacianFilter *New();
	vtkTypeMacro(vtkSignedDistanceFieldLaplacianFilter,vtkSimpleImageToImageFilter);

	// Description: 
	// The number of iterations
	vtkGetMacro(Iterations,int);
	vtkSetMacro(Iterations,int);

protected:

	vtkSignedDistanceFieldLaplacianFilter();
	~vtkSignedDistanceFieldLaplacianFilter() {};

	virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);

private:
	vtkSignedDistanceFieldLaplacianFilter(const vtkSignedDistanceFieldLaplacianFilter&);  // Not implemented.
	void operator=(const vtkSignedDistanceFieldLaplacianFilter&);  // Not implemented.

	//! Fill visit vector
	void CreateVisitOrder(int dim[3]);

	//! 
	void CreateLaplaceCoefficients();

	//! Fill coefficient cube with coeffiencts from a single Laplacian
	/** Translated by x,y,z and scaled by fac */
	void LocalLaplaceCoefficients(int x, int y, int z, double fac);

	//! Create filter coeffiencts lookup
	void CreateCoefficients();

	//! Simple visualisation
	void VisualiseCoefficients();

	struct CVoxelID
	{
		// voxel id (x,y,z)
		int id[3];
	};

	struct CFilterCoefficient
	{
		//! Coeffiecient displacement
		int p[3];

		//! Value
		double v;
	};

	//! Vector that list the visit order of the voxels
	std::vector<CVoxelID> m_VisitOrder;

	//! Do local filtering using a smoothing kernel
	void LocalSmoothingFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx);

	//! Do local filtering using a double Laplacian kernel
	void LocalDoubleLaplacianFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx);

	//! Cube used to calculate double Laplace coefficients
	std::vector<std::vector<std::vector<double> > > m_CoeffCube;

	//! A vector of the coeffiecients
	std::vector<CFilterCoefficient> m_FilterCoefficients;

	//! Number of iterations
	int Iterations;
};

#endif
