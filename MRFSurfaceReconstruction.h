#ifndef _MRFSurfaceReconstruction_h_
#define _MRFSurfaceReconstruction_h_

#include <string>
#include <iostream>

class vtkPoints;
class vtkPolyData;
class vtkImplicitFunction;
class vtkImageData;
class vtkImplicitVolume;
class C3dMDCamera;
class vtkCellLocator;

//! Encapsulates the functionality of the surface reconstruction algorithms in vtkImplicitLib
/** */
class CMRFSurfaceReconstruction
{
public:
	//! Default constructor
	CMRFSurfaceReconstruction();

	//! Destructor
	~CMRFSurfaceReconstruction();

	class CParms
	{
	public:

		CParms();

		//! Write parametes to file
		void DumpToFile(const std::string &name = "");

		//! Set input type
		/** Check the doc file for the normal path of processing for each type 
		0 : Points x,y,z                    | points with no normals 
		1 : Points with approximate normals | Normals can be estimated roughly by an acquisition device 
		2 : Points with good normals        | Points can be noise/outliers but still have good normal estimates 
		3 : Mesh                            | Mesh with connectivity but no normals 
		4 : Mesh with approximate normals   | Mesh with connectivity and rough estimates of normal directions 
		5 : Mesh with good normals          | Mesh with connectivity and well defined normals 
		*/
		void SetInputType(int type);

		//! Input file name
		std::string InputName;

		//! Short input name (without path and extension)
		std::string InputNameShort;

		//! Where is input placed
		std::string InputDir;

		//! Where to put output
		std::string OutPutDir;

		//! How much information should be written to disk
		/** 0: Nothing
		    1: Basics
			2: More than the basic..
			3: and so on... */
		int WriteLevel;

		//! Should the surface be computed and remeshed
		/** Sometimes it is only the distance field that is interesting and not the reconstructed surface.
		    0: No surface computed
			1: Surface computed and remeshed*/
		int DoComputeSurface;

		//! Place to find calibration files
		std::string CalibrationDir;

		//! The type of input data.
		/** Check the doc file for the normal path of processing for each type 
		0 : Points x,y,z                    | points with no normals 
		1 : Points with approximate normals | Normals can be estimated roughly by an acquisition device 
		2 : Points with good normals        | Points can be noise/outliers but still have good normal estimates 
		3 : Mesh                            | Mesh with connectivity but no normals 
		4 : Mesh with approximate normals   | Mesh with connectivity and rough estimates of normal directions 
		5 : Mesh with good normals          | Mesh with connectivity and well defined normals 
		*/
		int InputType;

		//! Estimate normals from input point using local PCA. Also removes noise and outliers 
		/** M1 switch from doc */
		bool EstimateNormals;

		//! Analyse the connected components of the point cloud
		/** M2 switch from doc. Can remove components and also create consistent normals directions */
		bool ConnectedComponentAnalysis;

		//! Use input data normals to create consistent normal directions
		bool UseReferenceNormalsInVoting;

		//! Scale input to a normalised size
		bool ScaleInput;

		//! Radius used when estimating normals
		double NormalRadius;

		//! Radius used when estimating normals (factor when using adaptive parms)
		double NormalRadiusFactor;

		//! 0 = search in the search radius for points, 1 = search for NormalPointsPerNormal
		int NormalSearchMode;

		//! Number of points to search for when estimating normals
		int NormalPointsPerNormal;

		//! Maximum value of 3rd eigenvalue when estimating normals (in percent)
		double Max3rdEigenval;

		// Threshold on the maximum distance from the point to the estimated plane
		double MaxPlaneDistance;

		//! Minimum size (number of points) of cliques that are accepted in the normal connectivity analysis
		int MinCliqueSize;

		//! When using adaptive parameters the MinCliqueSize is estimated as MinCliquePercent * TotalNumberOfPoins
		double MinCliquePercent;

		//! If set to true, only largest clique is kept, regardless of MinCliqueSize
		int LargestCliqueOnly;

		//! Distance used to determine if two points are neighbours when cliques are constructed
		double CliqueNeighbourDistance;

		//! Smallest cell size
		double MinCellSize;

		//! Set manual sample space. Default -1;
		double SampleSpace;

		//! The global factor used in calculating the alpha weight. Control the "spread"
		double kappa;

		//! When automatically computing sample spacing, this factor can be used to create finer resolution
		/** In the default behaviour this factor is multiplied the default minimum cells size. So:
		    SampleFactor < 1 : Smaller output cells
			SampleFactor > 1 : larger output cells.
			(default 2) */
		double SampleFactor;

		//! Pad voxels - how many voxels should be on each side of the point cloud in the volume
		int PadVoxels;

		//! Should the bounding volume be automatically reduced per iteration
		bool ShrinkBoundingVolume;

		//! Distance mode using when creating the distance map
		/** 0 = Nearest point Euclidean distance, 
		1 = Nearest point projected distance, 
		2 = Median of nearest points projected distance (L1 norm),
		3 = Average of nearest points projected distance (L2 norm)*/
		int DistanceMode;

		//! Is it a signed distance (1) or unsigned (0)
		/** The signed distance is the default. */
		int SignedDistance;

		//! How many local points should be used to calculate the median or the average
		int NumberOfDistances;
		
		//! Type of prior energy
		/** 0 = Difference of Laplacians, 1 = Difference of voxel values */
		int PriorType;

		//! Should local weighting be used
		bool UseLocalWeights;

		//! Max allowed iterations
		int MaxIterations;

		//! Global weight
		/** Beta = 0 : Pure prior (smoothing)
			Beta = 1 : Uses locals weight 100% */
		double GlobalBeta;

		//! Set if the ICM should only work in a band around the zero level after 2 levels
		/** Set to false if artefacts happens in the output */
		int BandedICM;

		//! Which polygoniser to use
		/** 0 - Bloomenthal
		    1 - vtkContourFilter (Marching Cubes)*/
		int Polygoniser;

		//! How many times smaller than the voxel spacing should the Bloomenthal cell be
		double CellSizeFactor;

		//! 0 = lookup in PD data (f.ex. local density), 1 = Euclidian squared distance to nearest PD point
		int WeightMode;

		//! Max distance used when computing the local weights using the distance-mode
		double LocalWeightMaxDist;

		//! If true, a patch will be cut out from the point cloud and the interpolation abilities of the algorithm is este
		bool UseLeavePatchOut;

		//! If true, noise is added to the point cloud (not compatible with UseLeavePatchOut)
		bool AddNoise;

		//! Nature of noise (0 = 5% outlier noise, 1 = Gaussian noise)
		int NoiseNature;

		//! Use multi level solving
		bool MultiLevel;

		//! The maximum number of voxels in the final volume in multi-level mode
		int MaxVolumeSize;

		//! Remesh to get better triangles
		/** 0: no remeshing, 1: classical remeshing, 2: only project point to zero level */
		int Remesh;

		//! Remesh to target edge lengths (edge lenghts specified as input data scalars)
		bool RemeshToTargetEdgeLengths;

		//! Clone textures from VRML to remeshed surface
		bool CloneTextures;

		//! Type of optimisation
		/** 0 : ICM
			1 : Conjugate Gradient
			2 : Sparse Cholesky */
		int Optimisation;

		//! If this set to true the following parameters are estimated based on the input point cloud statistics
		/** 
			NormalRadius
			MaxPlaneDistance
			MinCliqueSize
			LocalWeightMaxDist
			*/
		bool AdaptiveParameters;

		//! Compute the aggresively cropped surface
		bool AggresiveCrop;

		// Default iso-value: only used in unsigned distance field
		double IsoValue;
	};

	//! Compute in many steps and dump results to disc
	void Compute(CParms &parms);

	//! Compute the surface based on input data
	/** Return surface in output */
	bool OneShotCompute(CParms &parms, vtkPolyData *input, std::string &msg);

	//! Polygonise implicit volume and remesh it if possible
	/** Returns 0 in case of errors. 1 in case of ok, 2 in case of warning */
	static int PolygoniseAndRemesh(vtkImageData * SDF, vtkImplicitVolume * impvol, vtkPolyData * output, int Polygoniser, double CellSizeFactor, int remesh,
		vtkPolyData *targetEdgeLengthPD, bool FixNonManifold, double IsoValue, std::string &msg);

	//! Color a surface based on input data and a reconstructed surface
	static vtkPolyData *DoApproximateColorSurface(vtkPolyData *RawData, vtkPolyData *Surface, const CParms& parms);

	//! Utility function
	static vtkPolyData *CropSurfaceExternal(vtkPolyData *RawData, vtkPolyData *Surface, const CParms& parms);

	//! Compute implicit surface based on input data.
	/** Return implicit volume in output */
	bool OneShotComputeDistanceField(CParms &parms, vtkPolyData *input, std::string &msg );

	//! Compute normals from PCA (no surface reconstruction here)
	bool OneShotComputePCANormals(CParms &parms, vtkPolyData *input, std::string &msg );

	//! Data with estimated normals
	vtkPolyData *GetNormalData() const;

	//! The extracted surface
	vtkPolyData *GetMRFSurface() const;

	//! The cropped surface
	vtkPolyData *GetMRFSurfaceCropped() const;

	//! The aggresively cropped surface (extrapolation parts cropped away)
	vtkPolyData *GetMRFSurfaceAgressivelyCropped() const;
	
	//! The regularised distance field
	vtkImageData *GetMRFDistanceField() const;

	//! Get parameters that might have been set adaptively
	CParms GetUpdatedParms() const;

	//! Compute adaptive parameters
	bool AdaptiveParameters(vtkPolyData *pd, CParms &parms);

private:

	//! Read calibration files
	void ReadCalibrationFiles();

	//! Convert raw data
	bool ConvertRawScanData();

	//!
	bool CreateNormals();

	//! Visualise spheres used to cutting
	void CreateCutSpheres();

	//! Cut out data with sphere
	void CutWithSpheres();

	//! Add some noise
	void AddNoiseToCloud();

	//! Transform so pointcloud fills bounding box better
	void TransformToOptimalBoundingBox();

	//! Compute some statistiscs on input data
	/** Also set parameters if the adaptive parameter settings is chosen */
	bool InputDataStatistics();

	
	//! Create density map
	void CreateDensity();

	//!
//	void MRFPolygoniser();

	//! Crop the MRF surface
	bool CropSurface();

	//! Color surface based on the local weights
	/** The full weight volume is used as lookup */
	void WeightColorSurface();

	//! Approximate local weight by calculating the distance to each sample point
	void ApproximateWeightColorSurface(int mode);

	//! Compute the distance map based on the oriented point set
	void ComputeDistanceMap();

	//! MRF based regularisation of the distance map
	void MRFRegularisation();

	//! Computes the distance field and solves the MRF using a Sparse Cholesky solver
	void SparseCholeskyMRFRegularisation();

	//! Computes the distance field and solves the MRF in several levels
	bool MultiLevelMRFRegularisation();

	//! Estimate a new volume 
	/** Shrink it so it fits tightly */
	int EstimateNewVolumeDimensions(vtkImageData *vol, vtkPolyData *points, int *extent, double &sampleSpacing, double *origin, int maxSize, double MinCellSize);

	//! Estimate a new volume 
	/** Do not shrink it */
	int EstimateNewVolumeDimensionsWithoutShrink(vtkImageData *vol, vtkPolyData *points, int *extent, double &sampleSpacing, double *origin, int maxSize, double MinCellSize);

	// Check if edges of suggested new volume cuts iso-zero surface (there is a sign change at the edge)
	void CheckVolumeIfIsoSurfacePresentAtEdges(vtkImageData * vol, vtkPolyData * points, int * extent, double & sampleSpacing, double * origin, int maxSize, double MinCellSize);

	bool CheckYSliceIsoSurface(vtkImageData * vol, int y);

	bool CheckXSliceIsoSurface(vtkImageData * vol, int x);

	bool CheckZSliceIsoSurface(vtkImageData * vol, int z);

	// Get the bounds of a possible iso-surface in a volume. If the iso-surface extends beyond an edge that edge value is set to -1
	void GetIsoSurfaceBounds(vtkImageData * vol, double * bounds);
	
	//! Write a VTK file with the bounds of the image data
	void WriteBoundsFile(vtkImageData *dat, int iteration);

	//! Write diagnostic information
	void WriteDiagnosisFile(double SplineEnergy, int iteration);
	
	
	//! Write a VTK file with the voxel size
	void WriteVoxelSizeFile(vtkImageData *dat, int iteration);
		
		//! Refine an already computed MRF volume
	/** The Sparse Cholesky solver could be used to create the initial guess */
	void MultiLevelMRFRefinement();
	
	//! Extract iso surface from the original distance map
	void PolygoniseRawDistanceMap();

	//! Extract iso surface from the processed distance map
	bool PolygoniseDistanceMap();

	//! Calculate surface statistics
	bool SurfaceStatistics();

	//! Test of upscaling
	void PolygoniseUpscaledDistanceMap();
		
	//! Compute accuracy by comparing input points with the surface
	void ComputeAccuracy();

	//! Compute the accuracy of a cut out patch by comparing the points with the surface that interpolated
	void ComputePatchAccuracy();

	//! Check surface for consistency like existence of non-manifold edges
	/** returns true if surface if fine */
	bool CheckSurfaceConsistency();

	//! Remesh to getter better triangles
	bool Remesh();

	//! Clone texture from VRML file
	void CloneTexturesFromVRML();

	//! Use calibration files to recompute texture coordinate
	void ReTextureSurface();

	//! Check if a given point can be seen from a given camera
	bool CheckCameraVisibility( C3dMDCamera &cam, double * p, double *n, vtkCellLocator * locator, double &normScore );

	//! Create a polydata where alle vertices have been split
	vtkPolyData * SplitVertices(vtkPolyData *pd);

	//! Create a mesh, where each triangle only contains vertices seen from one camera
	void RelabelVisibilityMesh(vtkPolyData *pd);


	//! Calculate the center of mass of a mesh
	static void CenterOfMass(vtkPoints* pd, double *CM);

	//! Helper function to check if a file exists
	static bool file_exists(char const* fn);

	//! Delete all allocated objects
	void DeleteAll();

	//! Read texture images and create one combined texture
	bool CreateCombinedTexture(const std::string &textureName1,const std::string &textureName2, const std::string &TextureOut);

	//!
	CParms m_Parms;

	//! The raw input data
	vtkPolyData *RawData;

	//! Data with estimated normals
	vtkPolyData *NormalData;

	//! The distance field
	vtkImageData *MRFDistanceField;

	//! The extracted surface
	vtkPolyData *MRFSurface;

	//! The cropped surface
	vtkPolyData *MRFSurfaceCropped;

	//! The aggresively cropped surface
	vtkPolyData *MRFSurfaceAggresivelyCropped;
};

#endif
