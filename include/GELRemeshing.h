#ifndef _GELRemeshing_h_
#define _GELRemeshing_h_

#include <iostream>
#include <fstream>
#include <string>

#include <GEL/HMesh/Manifold.h>
#include <vector>

class vtkPolyData;
class vtkImplicitFunction;
class vtkImageData;
class vtkPointLocator;
class vtkDoubleArray;
class vtkCellLocator;

//! GEL based remeshing
/** The methods is an implementation of Mario Bosch and Leif Kobbelts Pliant remesher

	Implementation tips are from their SGP 2008 course notes.*/
class CGELRemeshing
	{
	public:
		//! Default constructor
		CGELRemeshing();

		//! Destructor
		virtual ~CGELRemeshing();

		bool Polygonise(vtkImageData * SDF, vtkImplicitFunction *func);

		bool Remesh(vtkImplicitFunction *func, vtkPolyData *pdin, vtkPolyData *pdout, double targetEdgeLength = -1, double edgeFactor = 1);

		// For unsigned distance fields. Here the closest point on the surface is found using bisection search along the normal
		bool RemeshWithBisectionProjection(vtkImplicitFunction* func, vtkPolyData* pdin, vtkPolyData* pdout, double targetEdgeLength = -1);

		// For unsigned distance fields. Many points are projected to the zero level. Mesh connectivity thrown away
		bool RemeshWithMultiPointProjection(vtkImplicitFunction* func, vtkPolyData* pdin, vtkPolyData* pdout, double targetEdgeLength = -1);
		
		bool RemeshDirectMesh(vtkPolyData *pdin, vtkPolyData *pdout, double targetLength);

		
		bool RemeshWithCurvature(vtkImplicitFunction *func, vtkPolyData *pdin, vtkPolyData *curv, vtkPolyData *pdout);

		// Remesh so each vertex length is close to a specified target length
		/* Works directly on the source mesh - meaning no implicit volume needed */
		bool RemeshWithTargetLengthsDirectlyOnMesh(vtkPolyData *pdin, vtkPolyData *curv, vtkPolyData *pdout);

		//! Only projecting to zero-level
		bool ProjectOnly(vtkImplicitFunction* func, vtkPolyData* pdin, vtkPolyData* pdout);


		//! Remeshing following Jakob Andreas Bærentzens implementation
		// bool RemeshJAB(vtkImplicitFunction *func, vtkPolyData *pdin, vtkPolyData *pdout);
		
		// Detect and visualise bad face normals (a bad face normal has an opposite normal to its neighbors)
		// Returns the number of badnormals
		static int VisualiseAdjacentFlippedNormals(vtkPolyData *pdin, vtkPolyData *pdout, double NormalLength);

		
		// Assume that input polydata is sampled on a sphere with a radius and centre 0,0,0
		bool RemeshSphere( vtkPolyData *pdin, vtkPolyData *pdout, double radius);

		//! Write some edge statistics to file
		static void DumpEdgeStatisticsToFile(vtkPolyData *pdin, const std::string& fname);
		
		//! Compute edge statistics
		static void ComputeEdgeStatistics(vtkPolyData *pdin, std::vector<double> &lengths);

		static double ComputeAverageEdgeLength(vtkPolyData *pdin);

		//! Dumpe min angle statistics to file
		static void DumpMinAngleStatisticsToFile( vtkPolyData *pdin, const std::string& fname );

		//! Write some triangle area statistics to file
		static void DumpAreaStatisticsToFile(vtkPolyData *pdin, const std::string& fname);

		// Should be moved to a general lib
		static bool VTK2GEL(vtkPolyData *pd, HMesh::Manifold &mesh);

		// Generate GEL mesh by adding faces one by one and stitching in the end. 
		// Should be more robust for non-manifold meshes
		static bool VTK2GELFaceByFace(vtkPolyData *pd, HMesh::Manifold &mesh);

		// Should be moved to a general lib
		static void GEL2VTK(HMesh::Manifold &mesh, vtkPolyData *pd);

		// SomeGELBoundaryTest
		bool TestGELBoundaryOps(vtkPolyData *pdin);

		bool TestGELBoundaryCollapseEdge( vtkPolyData *pdin );
		
	private:

		void SplitLongEdges(HMesh::Manifold &mesh, double high);

		//! A port of the GELVCECTOR version
		void SplitLongEdgesReengineered(HMesh::Manifold &mesh, double high);

		void CollapseShortEdges(HMesh::Manifold &mesh, double low, double high);
		
		void CollapseShortBoundaryEdges(HMesh::Manifold &mesh, double low, double high);

		void CollapseShortBoundaryEdgesWithCurvature(HMesh::Manifold &mesh, vtkPointLocator * locator, vtkDoubleArray *scalars);

		void CollapseShortEdgesWithCurvature(HMesh::Manifold &mesh, vtkPointLocator * locator, vtkDoubleArray *scalars);

		void SplitLongEdgesWithCurvature( HMesh::Manifold &mesh, vtkPointLocator * locator, vtkDoubleArray *scalars);

		
		void RemoveCaps(HMesh::Manifold &mesh);

		void EqualizeValences(HMesh::Manifold &mesh);

		void MaximiseMinAngle(HMesh::Manifold &mesh);
		
		
		void TangentialRelaxation(vtkImplicitFunction *func, HMesh::Manifold &mesh);

		void TangentialRelaxationOnlyMesh(HMesh::Manifold &mesh);
		
		//! Use the one ring area for each vertex as weight
		void AreaWeigthedTangentialRelaxation(vtkImplicitFunction *func, HMesh::Manifold &mesh);
		

		void AreaWeigthedTangentialRelaxationOnSphere(HMesh::Manifold &mesh);
		
		
		void MoveToVertexCenterMass(vtkImplicitFunction *func, HMesh::Manifold &mesh);
				

		void BoundaryRelaxation(vtkImplicitFunction *func, HMesh::Manifold &mesh);
		
		void BoundarySmoothing(HMesh::Manifold &mesh);
		
		void ProjectToSurface(vtkImplicitFunction *func, HMesh::Manifold &mesh, double SanityJump, int iterations = 1000);
		
		void ProjectToSurfaceByBisection(vtkImplicitFunction* func, HMesh::Manifold& mesh, double epsilon);


		void ProjectToSurfaceOnMesh(HMesh::Manifold &mesh, vtkCellLocator *cLocator, double SanityJump);

		void ProjectToSurfaceOnSphere( double radius, HMesh::Manifold &mesh, double SanityJump );

		// Visualise normals of a face - mostly debugging use for see if faces are flipped
		void VisualiseFaceNormals(HMesh::Manifold &m, std::string &oname);

		// Visualise gradients
		void VisualiseGradients(HMesh::Manifold& m, vtkImplicitFunction* func, std::string& oname);

		// Detect and visualise bad face normals (a bad face normal has an opposite normal to its neighbors)
		// Returns the number of badnormals
		int VisualiseBadFaceNormals(HMesh::Manifold &m, std::string &oname, bool writeFile);

		// Show the edges that has sharp angle 
		void VisualiseAngleBasedBoundary(HMesh::Manifold& m, std::string& oname, bool writeFile);


		//! The summed area surrounding a vertex
		double VertexArea(HMesh::Manifold &mesh, HMesh::VertexID vi);

		// Check if one edge surrounding a vertex is an edge. Both based on connectivity and sharp edges
		bool IsBoundaryWithAngleCheck(HMesh::Manifold& mesh, HMesh::VertexID vi);
//		void LaplacianSmooth(HMesh::Manifold& m, double weight, int max_iter);


		// Modified version of GEL TAL smoothing. Uses another boundary check and other things
		void TAL_smoothing_RAPA_Version(HMesh::Manifold& m, float w, int max_iter);

		bool ClosestPointStatistics(vtkPolyData* pd, double &avgDist);

};

#endif
