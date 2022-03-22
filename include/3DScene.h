#ifndef _3DScene_h_
#define _3DScene_h_

#include <string>
#include <vector>

#include "SurfaceProperties.h"
// #include "ColorManager.h"
#include "Process3DData.h"
#include "SumatraSettings.h"

class vtkActor2D;
class vtkPointLocator;
class vtkRenderWindowInteractor;
class vtkScalarBarActor;
class vtkAxesActor;
class vtkPlaneWidget;
class vtkSphereWidget;
class vtkTextActor;
class vtkPointPicker;
class vtkSphereSource;
class vtkPolyDataMapper;

//! Keep tracks of 3D objects and have entries to manipulation routines
class C3DScene
{
	public:
		//! Default constructor
		C3DScene(CSumatraSettings *Settings);

		//! Destructor
		virtual ~C3DScene();

		//! Initialise internal states
		void Init();

		//! Set if gradient background should be used
		void SetGradientBackground(bool flag);

		//! Get if gradient background is used
		bool GetGradientBackground() const;

		//! Remove all surfaces
		void RemoveAllSurfaces();

		//! Remove surface with id and move shuffle all down.
		void RemoveSurface(unsigned int objID);

		//! Get the renderer
		vtkRenderer *GetRenderer() const;

		//! Read a file with the given name and add it to the scene
		/** Returns false in case of errors */
		bool ReadFile(const std::string &fname);

		//! Save a file 
		/** Returns false in case of errors 
		    \param ApplyUserTransform Transform object first. */
		bool SaveFile(unsigned int id, const std::string &fname, bool ApplyUserTransform = false, bool Binary = false,
			bool WriteNormals = false, bool WriteScalars = false);

		//! Set if the axes should be shown
		// void SetShowAxes(bool flag);

		//! Reset camera
		void ResetCamera();

		//! Get number of surfaces
		int GetNumberOfSurfaces() const;

		//! Get name of a surface
		std::string GetSurfaceShortName(unsigned int id) const;

		//! Get name of a surface
		std::string GetSurfaceFullName(unsigned int id) const;

		//! Get if undo is available
		bool GetUndoAvaible(unsigned int id) const;

		//! Get actor from a surface
		vtkActor *GetActor(unsigned int id) const;

		//! Get mapper
		vtkMapper *GetMapper(unsigned int id) const;

		//! Set scalar range
		void SetScalarRange(unsigned int id, double *range);
	
		//! Remove scalars from object
		void RemoveScalars(unsigned int id);

		//! Remove normal data from object
		void RemoveNormals(unsigned int id);
		
		//! Returns false if no scalars
		/** Range is the range of the scalar values in poly data */
		bool GetSurfaceScalarRange(unsigned int id, double *range);

		//! Returns false if no scalars
		//! Range is the range in the mapper
		bool GetScalarRange(unsigned int id, double *range);

		//! Returns some computed surface values in string format
		std::string GetSurfaceValues(unsigned int id) const;

		//! Do ICP alignment of two surfaces
		/** \param DoICP if false the surfaces are compared without ICP alignment
		    \param ApplyTransformTo (0 = source surface, 1 = Copy of source, 2 = other surface) */
		//void DoICPAlignment(unsigned int SourceID, unsigned int TargetID, bool DoICP, 
		//	unsigned int ApplyTransformTo, unsigned int otherSurface, bool SignedDistance, bool MatchCentroids, int MaxLandmarks,
		//	bool CalculateRegErrors);

		//! Do ICP alignment
		bool DoICPAlignment(CICPParameters& parms);

		//! Do compare
		bool DoCompare(unsigned int SourceID, unsigned int TargetID, bool SignedDistance, bool ReplaceSource, bool ExcludeEdgePoints);
		
		//! Do project surface
		void DoProjectSurface(unsigned int SourceID, unsigned int TargetID, bool CreateDisplacement, bool ReplaceSource);
		
		
		//! Animate surface in targetID using a linear combination of SourceID and DisplacementID
		/** To start the animation call it with TargetID = 100000, a new id will be returned that will be used in the
		    following */
		void Animate(unsigned int SourceID, unsigned int DisplacementID, unsigned int &TargetID, double t);

		//! Do merge
		void DoMerge(unsigned int SourceID, unsigned int TargetID);

		//! Flip object (f.ex. changing right ear into left ear)
		void FlipObject(unsigned int SourceID, bool ReplaceSurface);

		//! Calculate normals for a surface
		/** \param ReplaceSource if true the source surface is replaced with the result */
		void CalculateNormals(unsigned int SourceID, bool ReplaceSurface, bool flip, bool split, double featureAngle);
		
		//! Calculate MRF Surface Reconstruction
		void CalculateMRFSurface(unsigned int SourceID, bool ReplaceSurface, bool RecomputeNorms, int priortype, 
			bool LargestCliqueOnly, bool MarchingCubes, bool ConComp, int InputType, int NumberOfPoints, int PPerNormal, int PPerDistance, double TriangleSizeFactor, bool UseTargetEdgeLengths);

		//! Remesh surface
		void DoRemesh(unsigned int SourceID, bool ReplaceSurface, bool ScalarTargetLengths, double UniformTargetLength);

		//! Calculate PCA normals
		void CalculatePCANormals(unsigned int SourceID, bool ReplaceSurface, bool ConComp, bool KeepLargestClique, bool UseOrgNormals, bool VisNormals, int AdaptiveOrManual, double RadiusFactor, double PointsPerNormal, double ManualNormalRadius, double ManualPlaneDist, int ManualMinCliqueSize, double CliqueNeigDist);

		//! Compute adaptive parameters for PCA normals
		void CalculatePCANormalsParameters(unsigned int SourceID, double RadiusFactor,  double &ManualNormalRadius, double &ManualPlaneDist, int &ManualMinCliqueSize, double &CliqueNeigDist);
				
		//! Compute the average edge lengths in a mesh
		bool ComputeAverageEdgelengths(unsigned int SourceID, double &avgL);

		//! Resample shape using random points. Loose connectivity
		void RandomResample(unsigned int SourceID, bool ReplaceSurface, int NumPoints, bool  CloneNormals);

		//! Subsample point cloud
		void SubSamplePointCloud(unsigned int SourceID, bool ReplaceSurface, int NumPoints);
		
		//! Removes points of the surface which is outside the interval given by thresholds
		/** \param ReplaceSource if true the source surface is replaced with the result */
		void ScalarThreshold(unsigned int SourceID, bool ReplaceSurface, double lowThreshold, double highThreshold);

		//! Smoothes the values of the scalars of the surface by replacing a scalar value with the average value of all scalars found in a sphere around the point
		/** \param ReplaceSource if true the source surface is replaced with the result */
		void ScalarFiltering(unsigned int SourceID, bool ReplaceSurface, double SphereRadius);


		//! Removes points of the surface which is outside the interval given by thresholds
		/** \param FindClosestPoints finds the closest point from source to each point in target within MaximumDistance from each point */
		void FindClosestPoints(unsigned int SourceID, unsigned int TargetID,bool ReplaceSurface, double MaximumDistance);

		//! Extract largest connected region
		/** \param ReplaceSource if true the source surface is replaced with the result 
		    \param RegionType (0 = all regions, 1 = largest region, 2 = outer region, 3 = closest point region, 4 = scalar connectivity)*/
		void Connectivity(unsigned int SourceID, bool ReplaceSurface, int RegionType, double *scalarRange, bool fullScalarMode,
			double *point);

		//! Subdivide surface
		/** \param subdivtype 0 = vtkLinearSubdivisionFilter, 1 = vtkButterflySubdivisionFilter,
		     2 = vtkLoopSubdivisionFilter */
		void SubdivideSurface(unsigned int SourceID, bool ReplaceSurface, int subdivtype, int subdivisions);

		//! Decimate surface
		/** \param decimtype 0 = vtkDecimatePro, 1 = vtkQuadricDecimation 
		    \param decimfactor Specify the desired reduction in the total number of polygons (e.g., if decimfactor is set to 0.9, this filter will try to reduce the data set to 10% of its original size).*/
		void DecimateSurface(unsigned int SourceID, bool ReplaceSurface, int decimtype, float decimfactor, bool preservetopology);

		//! Smooth surface
		/** \param smoothtype 0: Laplacian, 1: Sinc*/
		void SmoothSurface(unsigned int SourceID, bool ReplaceSurface, int smoothtype, int NumIt, double RelaxFactor,
			bool BoundarySmooth, bool FeatureEdgeSmooth, double FeatureAngle, bool GenerateErrScal);

		//! Create FastRBF surface
		/** */
		//void CreateFastRBFSurface(unsigned int SourceID, bool OutputSurface, bool OutputEstNormals,  bool OutputDensity,
		//	double MinNormalLength, double MaxNormalLength, double Accuracy, double Resolution,
		//	double EstNormalsRadius, double EstNormalsPlaneFactor, bool ErrorBarFitting, double ISOSmooth,
		//	bool DefineDensity, int InsideID, int OutsideID);

		//! Simply merge all sets into a point sets. Keeping normals and scalars
		void MergeAllSets();

		//! Close holes in surface
		void CloseHoles(unsigned int SourceID, bool ReplaceSurface);

		//! Does surface has normals
		bool SurfaceHasNormals(unsigned int id) const;

		//! Create a polygonal representation of surface normals
		/** Can also visualise normals where two adjacent faces have flipped normals. */
		void VisualiseNormals(unsigned int SourceID, bool VisualiseAdjacentFlipped);
		
		//! Create normals by using the camera position
		void NormalsFromCamera(unsigned int SourceID);
				
		//! Calculate feature edges
		/** \param EdgeType 0 = boundary, 1 = non-manifold, 2 = manifold, 3 = sharp */
		void ComputeFeatureEdges(unsigned int SourceID, bool boundary, bool nonManifold, bool manifold, bool sharp, double FeatureAngle);

		//! Set if the scalarbar is visible
		void SetScalarBarVisible(bool flag);

		//! Get if the scalarbar is visible
		bool GetScalarBarVisible() const;

		//! Set if the axes is visible
		void SetAxesVisible(bool flag);

		//! Get if the axes is visible
		bool GetAxesVisible() const;

		//! Initialise the plane widget
		void InitialisePlaneWidget(vtkRenderWindowInteractor *interactor, unsigned int objID);

		//! Initialise the sphere widget
		void InitialiseSphereWidget(vtkRenderWindowInteractor *interactor, unsigned int objID);

		//! Create a cylinder
		void CreateCylinder(double radius, double height, int resolution, bool capping, double* pos = NULL);

		//! Create a cube
		void CreateCube(double XLength, double YLength, double ZLength, double *pos = NULL);

		//! Create a trapezoid
		void CreateTrapezoid(double XLength, double YLength, double ZLength, double alpha, double beta, bool capping);

		//! Create a sphere
		void CreateSphere(double radius, double startTheta, double endTheta, double startPhi, double endPhi,
			int thetaRes, int phiRes, double *pos);
		
		//! Delete plane widget
		void DeletePlaneWidget();

		//! Delete sphere widget
		void DeleteSphereWidget();

		//! Cut given object with the plane
		void CutWithPlane(unsigned int objID);

		//! Create a mirror image of the selected object
		void MirrorWithPlane(unsigned int objID);

		//! Mark given object with sphere
		void MarkWithSphere(unsigned int objID);

		//! Create objects used in sphere marking
		void StartMarkWithSphere(unsigned int objID);

		//! Delete marker objects
		void EndMarkWithSphere(unsigned int objID);

		//! Flip the plane
		void FlipPlane();

		//! Undo the last operation done on the selected object
		/** This is currently only possible with cut planes */
		void UndoLastOperation(unsigned int objID);

		//! Save plane to file
		bool SavePlane(const std::string& fname) const;

		//! Load plane from file
		bool LoadPlane(const std::string& fname);

		//! Set the marker value
		void SetMarkerValue(double val);
		
		//! Get current marker value
		double GetMarkerValue() const;

		//! Set radius of sphere widget
		void SetSphereWidgetRadius(double val);

		//! Get radius of sphere widget
		double GetSphereWidgetRadius() const;
		
		//! Set marker text
		void SetStatusText(const std::string txt, bool visible);

		//! Rotate camera
		void RotateCamera();

		//! Pick a point and return description
		std::string PickPointAndGetText(double x, double y);

		//! Set the scalar lookup table to one of the predifined values
		void SetScalarLookupTableNum(int num);

		////! Set the default point size
		//void SetDefaultPointSize(int ps);

		//! Set the background color
		void SetBackgroundColor(int R, int G, int B);

		//! Get position of last picked point
		double* getLastPickedPoint();
	
	private:

		CSumatraSettings* mSettings = NULL;

		//! Initialise renderers etc
		void SetupScreenThings();

		//! Set up the scalar bar
		void SetupScalarBar();

		//! Update Lookup table
		// void UpdateLookupTable(vtkPolyData *pd);

		vtkRenderer *m_Renderer;

		//! A list of surfaces
		std::vector<CSurfaceProperties*> m_Surfaces;

		//! The scalar bar
		vtkScalarBarActor *m_scalarBar;

		//! Is the scalarbar visible
		bool m_ScalarBarVisible;

		//! The axes
		vtkAxesActor *m_Axes;

		//! Is the axes visible
		bool m_AxesVisible;

		//! Color lookup table
		vtkLookupTable *m_lookup;

		//! The plane widget used to cut
		vtkPlaneWidget *m_PlaneWidget;

		//! Sphere widget used to mark surfaces
		vtkSphereWidget *m_SphereWidget;

		//! The color manager used to get a new color for each shape
		// CColorManager colorManager;

		//! Add surface to renderer and do other minor things
		void AddSurfaceToRenderer(CSurfaceProperties *sp);

		//! Update all bounds
		void UpdateAllBounds();
	
		//! The bounds of all actor. (xmin,xmax, ymin,ymax, zmin,zmax) 
		double m_AllBounds[6];

		//! Locator used with marker
		vtkPointLocator *m_SphereMarkerLocator;

		//! Current marker value
		double m_MarkerValue;

		//! Actor for the status text
		vtkTextActor *m_StatusText;

		//! Actor for the gradient background
		vtkActor2D *m_BackgroundActor;

		//! Used for picking
		vtkPointPicker *m_PointPicker;

		//! To show where we picked
		vtkSphereSource *m_PickSphere;

		vtkPolyDataMapper *m_PickSphereMapper;

		vtkActor *m_PickSphereActor;

		//! Default pointsize should be changed into a general structure
		// int m_DefaultPointSize;

		//! Default ID of mirror surface
		int m_MirrorID;

		// Last picked point position
		double m_pickedPointPosition[3] = {0, 0, 0};

		C3DScene() {};
};

#endif
