#ifndef __vtkDistanceFieldMRFRegularisation_h
#define __vtkDistanceFieldMRFRegularisation_h

//#include "MRFEnergyComputer.h"
#include "vtkSimpleImageToImageFilter.h"
#include <vector>
//#include <vnl/vnl_cost_function.h>
//#include <vnl/vnl_vector.h>
//#include <vgl/vgl_point_3d.h>

class vtkIntArray;
class vtkDoubleArray;
class vtkPolyData;
class vtkIntArray;

#define VTK_WEIGHT_BASED_ON_INPUT 0
#define VTK_WEIGHT_BASED_ON_DISTANCE 1

#define VTK_PRIOR_ENERGY_DIFF_LAPLACE 0
#define VTK_PRIOR_ENERGY_DIFF_VOXEL   1
#define VTK_PRIOR_ENERGY_UNIFORMGRADIENT 2

#define VTK_OPTIMISATION_ICM 0
#define VTK_OPTIMISATION_CONJGRAD 1
#define VTK_OPTIMISATION_SPARSECHOL 2

// Performs regularisation of a signed distance field
// The approach is simalar to Markov Random Field optimisation
// Uses the original distance field in a regularisation scheme
class vtkDistanceFieldMRFRegularisation : public vtkSimpleImageToImageFilter
{
public:
	static vtkDistanceFieldMRFRegularisation *New();
	vtkTypeMacro(vtkDistanceFieldMRFRegularisation,vtkSimpleImageToImageFilter);

	// Description: 
	// The maximum number of iterations
	// It can still stop before if RMS is less than MinRMS
	vtkGetMacro(Iterations,int);
	vtkSetMacro(Iterations,int);

	// Description: 
	// Use ICM in a band
	vtkSetMacro(BandedICM, int);
	vtkGetMacro(BandedICM, int);
	vtkBooleanMacro(BandedICM, int);

	// Description: 
	vtkSetMacro(PerVoxelICMMemory, int);
	vtkGetMacro(PerVoxelICMMemory, int);
	vtkBooleanMacro(PerVoxelICMMemory, int);


	// Description: 
	// Global weight  alpha * originaldistance + (1-alpha) * smooth
	// alpha = 1: only original distance
	// alpha = 0: pure smoothing
	vtkGetMacro(GlobalBeta,double);
	vtkSetMacro(GlobalBeta,double);

	// Description: 
	// Stopping criteria. Stops when RMS less then MinRMS
	vtkGetMacro(MinRMS,double);
	vtkSetMacro(MinRMS,double);

	// Description: 
	// Conjugate gradient tolerance
	vtkGetMacro(CGTolerance,double);
	vtkSetMacro(CGTolerance,double);
	
	// Final energy using the spline term
	vtkGetMacro(FinalSplineEnergy,double);

	// Description: 
	//! The value used in normalising the distance weight
	vtkGetMacro(WeightValueHigh,double);
	vtkSetMacro(WeightValueHigh,double);

	// Description: 
	std::vector<double> GetValueTrack();

	// Set reference volume
	vtkSetMacro(ReferenceVolume,vtkImageData*);

	// Set reference polydata
	vtkGetMacro(ReferencePD,vtkPolyData*);
	vtkSetMacro(ReferencePD,vtkPolyData*);

	vtkSetClampMacro(WeightType,int,
		VTK_WEIGHT_BASED_ON_INPUT,VTK_WEIGHT_BASED_ON_DISTANCE);
	vtkGetMacro(WeightType,int);

	vtkSetClampMacro(PriorType,int,
		VTK_PRIOR_ENERGY_DIFF_LAPLACE,VTK_PRIOR_ENERGY_UNIFORMGRADIENT);
	vtkGetMacro(PriorType,int);

	// The type of optimisation used
	vtkSetClampMacro(Optimisation,int,
		VTK_OPTIMISATION_ICM,VTK_OPTIMISATION_SPARSECHOL);
	vtkGetMacro(Optimisation,int);

	// Get volume with weights
	vtkGetMacro(WeightVol,vtkImageData*);

	void SetWeightVol(vtkImageData* WV);

protected:

	vtkDistanceFieldMRFRegularisation();
	~vtkDistanceFieldMRFRegularisation();

	virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);

	void UseICMOptimisation( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input);
	
	void UseICMOptimisationMultithreaded( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input);
	
	void UseICMOptimisationWithSeperateOutputMultithreaded( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input);
	
	void UseICMOptimisationWithConvolution( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input);
	
	void UseICMOptimisationWithConvolutionMultithreaded( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input);
	
	

	

	

	double InnerLoopThread(int start, int stop, vtkDoubleArray * OutScalars, int * dims );
	
	struct ThreadArgs
	{
		int start;
		int stop;
		vtkDoubleArray *OutScalars;
		vtkDoubleArray *TempOut;
		int *dims;
	};

	double InnerLoopWithSeperateOutputThread(ThreadArgs Targs);

	double InnerLoopWithConvolutionThread(ThreadArgs Targs);

	double InnerLoopWithConvolutionBorderVersionThread(ThreadArgs Targs, bool topbottom);
		

	//! A memory volume is used to keep track of which voxels should be updated
	void UseICMOptimisationWithMemory( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input);
	
	
	void UseConjugateGradientOptimisation( int * dims, vtkDoubleArray * OutScalars);
	
	void UseSparseCholeskyOptimisation( int * dims, vtkDoubleArray * OutScalars );


private:
	vtkDistanceFieldMRFRegularisation(const vtkDistanceFieldMRFRegularisation&);  // Not implemented.
	void operator=(const vtkDistanceFieldMRFRegularisation&);  // Not implemented.

	//! Fill visit vector
	void CreateVisitOrder(int dim[3]);
	
	//! Fill visit vector. Fill it with a band around the zero set
	void CreateBandedVisitOrder(int dim[3], vtkImageData *infield);
		
	//! 
	void CreateLaplaceCoefficients();

	//! Create local weights base on reference volume and reference PD
	void CreateLocalWeights(vtkImageData* input);

	//! Weights based on distance to nearest point sample
	void DistanceBasedWeights( vtkImageData* input, vtkIntArray * RefVolData);
	
	//! Fill coefficient cube with coeffiencts from a single Laplacian
	/** Translated by x,y,z and scaled by fac */
	void LocalLaplaceCoefficients(int x, int y, int z, double fac);

	//! Create filter coeffiencts lookup
	void CreateCoefficients();

	//! Simple visualisation
	void VisualiseCoefficients();

	//! Helper function
	static double GetValue(vtkDoubleArray* data, int dims[3], int x, int y, int z);

	//! Helper function
	static void SetValue(vtkDoubleArray* data, int dims[3], int x, int y, int z, double val);

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
	double LocalSmoothingFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx);
	
	//! Do local filtering that tries to set the gradient size to 1
	double LocalUniformGradientFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx, double spacing);
	
	//! Do local filtering using a double Laplacian kernel
	/** Returns the change in local value */
	double LocalDoubleLaplacianFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx);
	
	//! Uses seperate output data
	double LocalDoubleLaplacianFilteringWithSeperateOutput(vtkDoubleArray* data, vtkDoubleArray* outdata, int dims[3], CVoxelID idx);

	//! Uses simple convolution
	double LocalDoubleLaplacianFilteringUsingConvolution(vtkDoubleArray* data, vtkDoubleArray* outdata, int dims[3],  int xo, int yo, int zo);
	
	//! Uses simple convolution and wraps border of image
	double LocalDoubleLaplacianFilteringUsingConvolutionAndBorderCheck(vtkDoubleArray* data, vtkDoubleArray* outdata, int dims[3], int xo, int yo, int zo );
			
	//! Do local filtering using a double Laplacian kernel. Banded version without boundary check
	/** Returns the change in local value */
	double LocalBandedDoubleLaplacianFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx);
	
	//! Use last change information to predict future
	void ForwardPrediction(vtkDoubleArray* data, vtkDoubleArray* changes, int dim[3]);

	//! Util function to calculate median
	void Median(std::vector<double>& x, double fractile, double& median);

	void CreateLookup(int * dims, std::vector<int> &Lookup, std::vector<vgl_point_3d<int> >& PosLU);

	//! A sparse matrix entry
	class CMatrixEntry
	{
	public:
		int i;
		int j;
		double v;
	};

	void CreateSparseMatrix( int * dims, std::vector<CMatrixEntry> &matrix, std::vector<double>& b );

	void FindMatrixEntriesForNeighbours(int * dims, std::vector<int> &Lookup, std::vector<vgl_point_3d<int> >& PosLU, int x, int y, int z, int offset, std::vector<CMatrixEntry>& entries, std::vector<double> &b);

	void WriteMatrixAndB(std::vector<CMatrixEntry>& matrix, std::vector<double>& b);

	void CheckMatrixSymmetry(std::vector<CMatrixEntry>& matrix, std::vector<double>& b);

	void ReadResult(vtkDoubleArray * OutScalars);

	//! Cube used to calculate double Laplace coefficients
	std::vector<std::vector<std::vector<double> > > m_CoeffCube;

	//! A vector of the coeffiecients
	std::vector<CFilterCoefficient> m_FilterCoefficients;

	//! The filter coefficients as a double vector - for fast access
	double *m_ConvolutionCoefficients;


	//! Number of iterations
	int Iterations;

	//! Original scalars
	vtkDoubleArray *OrgScalars;

	//! Global weight  alpha * originaldistance + (1-alpha) * smooth
	double GlobalBeta;

	//! See vtkOrientedPointSetDistanceFilter for a description
	vtkImageData* ReferenceVolume;

	//! The values in ReferenceVolume are indices into the ReferencePD
	/** This can be used to create weigthing functions */
	vtkPolyData* ReferencePD;

	//! Local alphas based on reference volume and reference polydata
	vtkDoubleArray *LocalWeights;

	//! Bookkeeping volume
	vtkIntArray *ICMUpdateVolume;

	//! Volume with local weights
	vtkImageData *WeightVol;

	bool DeleteWeightVol;

	//! Should local alphas or global alpha be used
	bool UseLocalWeighting;

	//! Stopping criteria
	double MinRMS;

	//! What type of weight should be used
	int WeightType;

	//! Type of prior energy
	/** 0 = Difference of Laplacians, 1 = Difference of voxel values */
	int PriorType;

	//! Track values in one voxel
	std::vector<double> ValueTrack;

	//! Check how many visits each voxel got
//	std::vector<int> DebugVolume;

	//! Keep track of voxels changes
	/** Can be used to forward predictions */
//	vtkImageData *ChangeVolume;

	//! How long should the prediction be
	double PredictionStep;

	//! The value used in normalising the distance weight
	double WeightValueHigh;

	//! Tolerance for conjugate gradient
	double CGTolerance;

	//! The type of optimisation
	int Optimisation;

	//! Use ICM in a band
	int BandedICM;

	//! Use the result of last iteration to determine which voxels should be updated
	int PerVoxelICMMemory;

	// How many threads can be used. -1 single threaded. 0 automatically determine number of threads (#cores-1)
	int NumThreads;

	// Final spline energy
	double FinalSplineEnergy;
};
//
////! Minimiser for MRF
//class CMRFCostFunction : public vnl_cost_function
//{
//public:
//	CMRFCostFunction(int dims[3], int nunknowns, vtkDoubleArray *OS, vtkDoubleArray *LW, double GB, int EnergyType);
//
//	void SetDimensions(int dims[3]);
//
//	double f( vnl_vector< double > const& x );
//
//	void gradf (vnl_vector< double > const &x, vnl_vector< double > &gradient );
//
//
//private:
//
//	CMRFEnergyComputer m_EnergyComputer;
//
//	// 0: Spline
//	// 1: Membrane
//	int m_EnergyType;
//
//	CMRFCostFunction();
//};

#endif
