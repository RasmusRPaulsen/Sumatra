#ifndef _Process3DData_h_
#define _Process3DData_h_

#include <string>
#include <vtkPolyData.h>
#include <vtkMatrix4x4.h>


//! Parameters used by the ICP transform
class CICPParameters
{
public:
	//! Constructor
	CICPParameters();

	//! ID of source surface
	unsigned int SourceID;

	//! ID of target surface
	unsigned int TargetID;

	//! ID of surface that the transform can be applied to
	unsigned int transSurfaceID;

	//! Which surface should the transform be applied to
	enum eICPApplyTransformTo {
		eSource,				//!< Apply transform to source (replace source surface)
		eSourceCopy,			//!< Apply transform to copy of source surface
		eOtherSurface			//!< Apply transform to other surface
	};

	//! Which surface to apply transform to
	eICPApplyTransformTo ApplyTransformTo;

	//! Should the registration errors be calculated as signed distances
	bool SignedDistance;

	//! Should centroids of the surfaces be matched initially
	bool MatchCentroids; 

	//! How many landmarks to use
	int MaxLandmarks;

	//! Maximum number of iterations
	int MaxIterations;

	//! Mean distance stopping criteria
	double MaxMeanDistance;

	//! Should point landing on the edge of the target be excluded (constrained ICP)
	bool ExcludeEdgePoints;

	//! Should the registration errors be calculated
	/** If yes, the scalar values of the transformed surface will be overwritten */
	bool CalculateRegErrors;
	
	//! Source surface
	vtkPolyData *source;

	//! Target surface
	vtkPolyData *target;

	//! Surface to be moved
	vtkPolyData *move;

	//! Transform of the source surface
	vtkMatrix4x4 *sourcePreTransform;

	//! Transform of the target surface
	vtkMatrix4x4 *targetPreTransform;

	//! Transform of the surface to be moved
	vtkMatrix4x4 *moveSurfacePreTransform;

	//! The transformed surface
	/** If CalculateRegErrors == true, the scalar field will contain registration errors after the ICP */
	vtkPolyData *transformedSurface;

//	vtkPolyData *surfaceErrors;
};

//! Process 3D data 
/** Routines to do some 3D Processing */
class CProcess3DData
{
	public:
		//! Default constructor
		CProcess3DData();

		//! Destructor
		virtual ~CProcess3DData();

		//! ICP align source with target surface
		/** A custom transform matrix can be given for both the source and the target surface.
		    These can for example come from an actor transform.
			The transformed source and a registration error surface is returned.
			\param DoICP if false the surfaces are compared without ICP alignment */
		static void DoICP(CICPParameters& parms);

		//! Subdivide surface
		/** \param subdivtype 0 = vtkLinearSubdivisionFilter, 1 = vtkButterflySubdivisionFilter,
		     2 = vtkLoopSubdivisionFilter */
		static void DoSubdivideSurface(vtkPolyData *source, int subdivisions, vtkPolyData *subdivided, 
			int subdivtype);

		//! Decimate surface
		/** \param decimtype 0 = vtkDecimatePro, 1 = vtkQuadricDecimation 
		    \param decimfactor Specify the desired reduction in the total number of polygons (e.g., if decimfactor is set to 0.9, this filter will try to reduce the data set to 10% of its original size).*/
		static void DoDecimateSurface(vtkPolyData *source, vtkPolyData *decimated, 
			int decimtype, float decimfactor, bool preservetopology);

		//! Smooth surface
		/** */
		static void DoSmoothSurface(vtkPolyData *source, vtkPolyData *smooth,int smoothtype, int NumIt, double RelaxFactor,
			bool BoundarySmooth, bool FeatureEdgeSmooth, double FeatureAngle, bool GenerateErrScal);

		//! Calculate feature edges
		/** \param EdgeType 0 = boundary, 1 = non-manifold, 2 = manifold, 3 = sharp */
		static void DoComputeFeatureEdges(vtkPolyData *source, vtkPolyData *FeatureEdges, int EdgeType, double FeatureAngle);

		//! Extract the outer surface from a polydata containing a set of distinct part
		/** This is a special case of connectivity filtration */
		static void ExtractOuterSurface(vtkPolyData *source, vtkPolyData *outsurf);

		private:
};

#endif
