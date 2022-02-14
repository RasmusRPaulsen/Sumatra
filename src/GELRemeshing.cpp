#include "GELRemeshing.h"
#include "GeneralUtils.h"
#include <vtkpolydata.h>
#include <vtkCellArray.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/cleanup.h>
#include <vtkCellLocator.h>
//#include <HMesh/HalfEdgeHandle.h>
//#include <HMesh/FaceCirculator.h>
//#include <HMesh/VertexCirculator.h>
//#include <HMesh/Vector.h>
//#include <HMesh/caps_and_needles.h>
#include <GEL/HMesh/refine_edges.h>
//#include <HMesh/caps_and_needles.h>
#include <GEL/HMesh/mesh_optimization.h>
#include <GEL/HMesh/triangulate.h>
#include <GEL/HMesh/smooth.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/Geometry/QEM.h>
#include <vtkPointLocator.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>


#include <vtkMath.h>
#include <vtkImplicitFunction.h>
#include <ImplicitFunctionUtils.h>
#include "vtkExtMisc.h"
#include <vtkImageData.h>
#include "BloomenthalPolygonizer.h"
#include "MeshMeasures.h"

#include <string>
#include <vector>
#include <GEL/CGLA/Vec3d.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include "MultiTimer.h"


CGELRemeshing::CGELRemeshing()
{
}

CGELRemeshing::~CGELRemeshing()
{
}

/*
The split long edges(high) function visits all edges of the current mesh. If an edge is longer
than the given threshold high, the edge is split at its midpoint and the two adjacent triangles
are bisected (2-4 split).

split long edges(high)
while exists edge e with length(e)>high do
split e at midpoint(e)
*/
void CGELRemeshing::SplitLongEdges( HMesh::Manifold &mesh, double high )
{
	std::cout << "Split";
	int splits = refine_edges(mesh, high);
	std::cout << " " << splits << ". ";
}

void CGELRemeshing::SplitLongEdgesReengineered( HMesh::Manifold &mesh, double high )
{
	std::cout << "Split";
//	std::cout << "Split with high " << high;
	int work = 0;

	std::vector<HMesh::HalfEdgeID> hedges;
	hedges.reserve(mesh.no_halfedges());
	std::copy(mesh.halfedges_begin(), mesh.halfedges_end(), std::back_inserter(hedges));

	HMesh::HalfEdgeAttributeVector<int> touched(mesh.allocated_vertices(), 0);

//	for(std::vector<HMesh::HalfEdgeID>::iterator h = hedges.begin(); h != hedges.end(); ++h)
	for (HMesh::HalfEdgeID h : hedges)
	{
		HMesh::Walker w = mesh.walker(h);
//		HMesh::Walker w = mesh.walker(*h);

		if(!mesh.in_use(h) || w.face() == HMesh::InvalidFaceID || length(mesh, h) < high || touched[h])
		{
			continue;
		}

		++work;

		touched[w.opp().halfedge()] = 1;
		HMesh::VertexID v = mesh.split_edge(h);

		HMesh::Walker wv = mesh.walker(v);
//		HMesh::Walker wv = mesh.walker(v);

		bool useTriangulation = false;

		if (useTriangulation)
		{
			// 28/1-2020 : This is probably faster - but is not tested enough
			HMesh::FaceID f1 = wv.opp().face();
			if (f1 != HMesh::InvalidFaceID)
				HMesh::triangulate(mesh, f1);

			HMesh::FaceID f2 = wv.prev().face();
			if (f2 != HMesh::InvalidFaceID)
				HMesh::triangulate(mesh, f2);
		}
		else
		{
			HMesh::FaceID f1 = wv.opp().face();
			if(f1 != HMesh::InvalidFaceID)
			   mesh.split_face_by_vertex(f1);

			HMesh::FaceID f2 =  wv.prev().face();
			if(f2 != HMesh::InvalidFaceID)
			    mesh.split_face_by_vertex(f2);
		}
	}
	std::cout << "(" << work <<  ")";
}


int CGELRemeshing::VisualiseAdjacentFlippedNormals(vtkPolyData *pdin, vtkPolyData *pdout, double NormalLength)
{
	int badNormals = 0;

	// Hard coded length - todo
	double NL = 1;
	if (NormalLength != 0)
		NL = NormalLength;

	HMesh::Manifold m;

	if (!VTK2GELFaceByFace(pdin, m))
		return false;

	vtkPoints *points = vtkPoints::New();
	vtkCellArray *lines = vtkCellArray::New();

	for (auto h : m.halfedges())
	{
		HMesh::Walker w = m.walker(h);
		HMesh::FaceID f1 = w.face();
		if (f1 == HMesh::InvalidFaceID)
		{
			//		std::cout << "Bad face ID" << std::endl;
			continue;
		}
		HMesh::FaceID f2 = w.opp().face();
		if (f2 == HMesh::InvalidFaceID)
		{
			continue;
		}

		CGLA::Vec3d n1 = HMesh::normal(m, f1);
		CGLA::Vec3d n2 = HMesh::normal(m, f2);

		// Check if normals points the same way
		if ((n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) > -0.5)
		{
			continue;
		}
		//		std::cout << "Opposing normals found between face " << f1 << " and " << f2 << std::endl;
		badNormals++;

		//		double p[3];
		double n[3];
		double pt[3];

		n[0] = n1[0];
		n[1] = n1[1];
		n[2] = n1[2];

		CGLA::Vec3d p(0);
		HMesh::Walker w2 = m.walker(f1);
		for (; !w2.full_circle(); w2 = w2.next())
		{
			p += m.pos(w2.vertex());
		}
		p /= w2.no_steps();
		//CGLA::Vec3d p1 = m.pos(w.vertex());
		//CGLA::Vec3d p2 = m.pos(w.next().vertex());
		//CGLA::Vec3d p3 = m.pos(w.prev().vertex());

		//ptransSource->GetOutput()->GetPoint(i, p);
		//normals->GetTuple(i, n);

		vtkMath::Normalize(n);
		n[0] *= NL;
		n[1] *= NL;
		n[2] *= NL;

		pt[0] = n[0] + p[0];
		pt[1] = n[1] + p[1];
		pt[2] = n[2] + p[2];

		double pt2[3];
		pt2[0] = p[0];
		pt2[1] = p[1];
		pt2[2] = p[2];

		lines->InsertNextCell(2);
		vtkIdType id = points->InsertNextPoint(pt2);
		lines->InsertCellPoint(id);
		id = points->InsertNextPoint(pt);
		lines->InsertCellPoint(id);
	}

	pdout->SetPoints(points);
	points->Delete();
	pdout->SetLines(lines);
	lines->Delete();

	return badNormals;
}

void CGELRemeshing::VisualiseAngleBasedBoundary(HMesh::Manifold& m, std::string& oname, bool writeFile)
{
	vtkPolyData* pd = vtkPolyData::New();
	vtkPoints* points = vtkPoints::New();
	vtkCellArray* lines = vtkCellArray::New();

	double SharpAngle = 180.0 - 30.0;
	//double athres = cos(SharpAngle * vtkMath::Pi() / 180.0);
	//double athres = cos(SharpAngle * vtkMath::Pi() / 180.0);
	double athres = -0.5;
	// std::cout << "Angle threshold: " << athres << std::endl;


	for (auto h : m.halfedges())
	{
		HMesh::Walker w = m.walker(h);
		HMesh::FaceID f1 = w.face();
		if (f1 == HMesh::InvalidFaceID)
		{
			//		std::cout << "Bad face ID" << std::endl;
			continue;
		}
		HMesh::FaceID f2 = w.opp().face();
		if (f2 == HMesh::InvalidFaceID)
		{
			continue;
		}

		CGLA::Vec3d n1 = HMesh::normal(m, f1);
		CGLA::Vec3d n2 = HMesh::normal(m, f2);

		// Check if normals points the same way
		if ((n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) > athres)
		{
			continue;
		}

		//double pt[3];
		//double pt2[3];
		CGLA::Vec3d p1 = m.pos(w.vertex());
		CGLA::Vec3d p2 = m.pos(w.prev().vertex());
		

		lines->InsertNextCell(2);
		vtkIdType id = points->InsertNextPoint(p1[0], p1[1], p1[2]);
		lines->InsertCellPoint(id);
		id = points->InsertNextPoint(p2[0], p2[1], p2[2]);
		lines->InsertCellPoint(id);
	}

	pd->SetPoints(points);
	points->Delete();
	pd->SetLines(lines);
	lines->Delete();

	if (writeFile)
	{
		vtkExtMisc::WritePDVTK(pd, oname);
	}

	pd->Delete();
}


int CGELRemeshing::VisualiseBadFaceNormals(HMesh::Manifold &m, std::string &oname, bool writeFile)
{
	int badNormals = 0;

	// Hard coded length - todo
	double NL = 1;

	vtkPolyData *pd = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *lines = vtkCellArray::New();

	for (auto h : m.halfedges())
	{
		HMesh::Walker w = m.walker(h);
		HMesh::FaceID f1 = w.face();
		if (f1 == HMesh::InvalidFaceID)
		{
	//		std::cout << "Bad face ID" << std::endl;
			continue;
		}
		HMesh::FaceID f2 = w.opp().face();
		if (f2 == HMesh::InvalidFaceID)
		{
			continue;
		}

		CGLA::Vec3d n1 = HMesh::normal(m, f1);
		CGLA::Vec3d n2 = HMesh::normal(m, f2);

		// Check if normals points the same way
		if ((n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) > -0.5)
		{
			continue;
		}
//		std::cout << "Opposing normals found between face " << f1 << " and " << f2 << std::endl;
		badNormals++;

		//		double p[3];
		double n[3];
		double pt[3];

		n[0] = n1[0];
		n[1] = n1[1];
		n[2] = n1[2];

		CGLA::Vec3d p(0);
		HMesh::Walker w2 = m.walker(f1);
		for (; !w2.full_circle(); w2 = w2.next())
		{
			p += m.pos(w2.vertex());
		}
		p /= w2.no_steps();
		//CGLA::Vec3d p1 = m.pos(w.vertex());
		//CGLA::Vec3d p2 = m.pos(w.next().vertex());
		//CGLA::Vec3d p3 = m.pos(w.prev().vertex());

		//ptransSource->GetOutput()->GetPoint(i, p);
		//normals->GetTuple(i, n);

		vtkMath::Normalize(n);
		n[0] *= NL;
		n[1] *= NL;
		n[2] *= NL;

		pt[0] = n[0] + p[0];
		pt[1] = n[1] + p[1];
		pt[2] = n[2] + p[2];

		double pt2[3];
		pt2[0] = p[0];
		pt2[1] = p[1];
		pt2[2] = p[2];

		lines->InsertNextCell(2);
		vtkIdType id = points->InsertNextPoint(pt2);
		lines->InsertCellPoint(id);
		id = points->InsertNextPoint(pt);
		lines->InsertCellPoint(id);
	}

	pd->SetPoints(points);
	points->Delete();
	pd->SetLines(lines);
	lines->Delete();

	if (writeFile && badNormals > 0)
	{
		vtkExtMisc::WritePDVTK(pd, oname);
	}

	pd->Delete();
	//
	//for (auto h : m.halfedges())
	//	{
	//		Walker w = m.walker(h);
	//		HMesh::FaceId f1 = w.face();
	//		FacedId f2 = w.opp().face();
	//	if (f2 != HMesh::InvalidFaceID)
	//		vec3d n1 = normal(m, f1);
	//		vec3d n2 = normal(m, f2);

	//	}
	return badNormals;
}


void CGELRemeshing::VisualiseGradients(HMesh::Manifold& m, vtkImplicitFunction* func, std::string& oname)
{
	//double NL = 1;

	vtkPolyData* pd = vtkPolyData::New();
	vtkPoints* points = vtkPoints::New();
	vtkCellArray* lines = vtkCellArray::New();

	for (HMesh::VertexID v : m.vertices())
	{
		CGLA::Vec3d p = m.pos(v);

		double pt[3];
		pt[0] = p[0];
		pt[1] = p[1];
		pt[2] = p[2];

		double* n = func->FunctionGradient(pt);
		double NL = func->FunctionValue(pt);

		vtkMath::Normalize(n);
		n[0] *= NL;
		n[1] *= NL;
		n[2] *= NL;

		double pt2[3];
		pt2[0] = n[0] + p[0];
		pt2[1] = n[1] + p[1];
		pt2[2] = n[2] + p[2];

		lines->InsertNextCell(2);
		vtkIdType id = points->InsertNextPoint(pt2);
		lines->InsertCellPoint(id);
		id = points->InsertNextPoint(pt);
		lines->InsertCellPoint(id);
	}

	pd->SetPoints(points);
	points->Delete();
	pd->SetLines(lines);
	lines->Delete();

	vtkExtMisc::WritePDVTK(pd, oname);

	pd->Delete();
}


void CGELRemeshing::VisualiseFaceNormals(HMesh::Manifold &m, std::string &oname)
{

	//// Get an estimate of the size of the pointcloud
	//double bounds[6];
	//ptransSource->GetOutput()->GetBounds(bounds);

	//double l = sqrt((bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
	//	(bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
	//	(bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));

	// Normal length. 1 % of bounding box diagonal
//	double NL = l * 0.01;

	// Hard coded length - todo
	double NL = 1;

	vtkPolyData *pd = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *lines = vtkCellArray::New();


	for (auto h : m.halfedges())
	{
		HMesh::Walker w = m.walker(h);
		HMesh::FaceID f1 = w.face();
		if (f1 == HMesh::InvalidFaceID)
		{
//			std::cout << "Bad face ID" << std::endl;
			continue;
		}

		CGLA::Vec3d n1 = HMesh::normal(m, f1);

//		double p[3];
		double n[3];
		double pt[3];

		n[0] = n1[0];
		n[1] = n1[1];
		n[2] = n1[2];

		CGLA::Vec3d p(0);
		HMesh::Walker w2 = m.walker(f1);
		for (; !w2.full_circle(); w2 = w2.next())
		{
			p += m.pos(w2.vertex());
		}
		p /= w2.no_steps();
		//CGLA::Vec3d p1 = m.pos(w.vertex());
		//CGLA::Vec3d p2 = m.pos(w.next().vertex());
		//CGLA::Vec3d p3 = m.pos(w.prev().vertex());

		//ptransSource->GetOutput()->GetPoint(i, p);
		//normals->GetTuple(i, n);

		vtkMath::Normalize(n);
		n[0] *= NL;
		n[1] *= NL;
		n[2] *= NL;

		pt[0] = n[0] + p[0];
		pt[1] = n[1] + p[1];
		pt[2] = n[2] + p[2];

		double pt2[3];
		pt2[0] = p[0];
		pt2[1] = p[1];
		pt2[2] = p[2];

		lines->InsertNextCell(2);
		vtkIdType id = points->InsertNextPoint(pt2);
		lines->InsertCellPoint(id);
		id = points->InsertNextPoint(pt);
		lines->InsertCellPoint(id);
	}

	pd->SetPoints(points);
	points->Delete();
	pd->SetLines(lines);
	lines->Delete();


	vtkExtMisc::WritePDVTK(pd, oname);

	pd->Delete();
	//
	//for (auto h : m.halfedges())
	//	{
	//		Walker w = m.walker(h);
	//		HMesh::FaceId f1 = w.face();
	//		FacedId f2 = w.opp().face();
//	if (f2 != HMesh::InvalidFaceID)
	//		vec3d n1 = normal(m, f1);
	//		vec3d n2 = normal(m, f2);

	//	}

}

void CGELRemeshing::SplitLongEdgesWithCurvature( HMesh::Manifold &mesh, vtkPointLocator * locator, vtkDoubleArray *scalars)
{
	std::cout << "Split";
	int work = 0;

	std::vector<HMesh::HalfEdgeID> hedges;
	hedges.reserve(mesh.no_halfedges());
	std::copy(mesh.halfedges_begin(), mesh.halfedges_end(), std::back_inserter(hedges));

	HMesh::HalfEdgeAttributeVector<int> touched(mesh.allocated_halfedges(), 0);

//	for(std::vector<HMesh::HalfEdgeID>::iterator h = hedges.begin(); h != hedges.end(); ++h)
	for (HMesh::HalfEdgeID h : hedges)
	{
		HMesh::Walker w = mesh.walker(h);

//		if(!mesh.in_use(*h) || w.face() == HMesh::InvalidFaceID || length(mesh, *h) < high || touched[*h])
		if(!mesh.in_use(h) || w.face() == HMesh::InvalidFaceID || touched[h])
		{
			continue;
		}
		// check target length
		double len = length(mesh, h);

		CGLA::Vec3d mid_point = 0.5 * (mesh.pos(w.vertex())+mesh.pos(w.opp().vertex()));

		double p[3];
		p[0] = mid_point[0];
		p[1] = mid_point[1];
		p[2] = mid_point[2];

		int id = locator->FindClosestPoint(p);
		double trgVal = scalars->GetValue(id);
		if (trgVal == 0)
		{
			std::cerr << "Target value of 0 found on target mesh - can not remesh to 0 length. Setting it to 1!!!!!!" << std::endl;
			trgVal = 1;
		}

		double high = 4.0/3.0 * trgVal;

		if (len < high)
		{
			continue;
		}

		++work;

		touched[w.opp().halfedge()] = 1;
		HMesh::VertexID v = mesh.split_edge(h);

		HMesh::Walker wv = mesh.walker(v);

		HMesh::FaceID f1 = wv.opp().face();
		if(f1 != HMesh::InvalidFaceID)
			mesh.split_face_by_vertex(f1);

		HMesh::FaceID f2 =  wv.prev().face();
		if(f2 != HMesh::InvalidFaceID)
			mesh.split_face_by_vertex(f2);

	}
	std::cout << "(" << work <<  ")";
}


double vertex_error(const HMesh::Manifold& m, HMesh::VertexID v, const CGLA::Vec3d& pb)
{
	Geometry::QEM q;
	CGLA::Vec3d pa(m.pos(v));

	for(HMesh::Walker vj = m.walker(v); !vj.full_circle(); vj = vj.circulate_vertex_cw()){
		HMesh::FaceID f = vj.face();
		if(f != HMesh::InvalidFaceID)
			q += Geometry::QEM(CGLA::Vec3d(0), CGLA::Vec3d(normal(m, f)));
	}
	return q.error(pb - pa);
}


/*
The collapse short edges(low, high) function collapses and thus removes all edges that are
shorter than a threshold low. Here one has to take care of a subtle problem: by collapsing
along chains of short edges the algorithm may create new edges that are arbitrarily long and
thus undo the work that was done in split long edges(high). This issue is resolved by testing
before each collapse whether the collapse would produce an edge that is longer than high. If so,
the collapse is not executed.

collapse short edges(low, high)
finished = false
 while exists edge e with length(e)<low and not finished do
 finished = true
  let e=(a,b) and let a[1],...,a[n] be the 1-ring of a
  collapse ok = true
  for i = 1 to n do
  if length(b,a[i])>high then
  collapse ok = false
  if collapse ok then
  collapse a into b along e
  finished = false
   */
void CGELRemeshing::CollapseShortEdges(HMesh::Manifold &mesh, double low, double high )
{
	std::cout << "Collapse";
	// We do not do it exactly as stated as above. We do not check if we create long edges
//	remove_needles(mesh, low);

	bool did_work = false;
	int collapses = 0;
	int PrevCol = 0;
	// remove needles until no more can be removed
	do
	{
		did_work = false;
		for(HMesh::VertexIDIterator v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
		{
			// don't attempt to remove needles if vertex is boundary
			// if(boundary(mesh, *v))
			if(IsBoundaryWithAngleCheck(mesh, *v))
			{
				continue;
			}

			for(HMesh::Walker vj = mesh.walker(*v); !vj.full_circle(); vj = vj.circulate_vertex_cw())
			{
				// don't attempt to remove needles if vertex of jumper halfedge is boundary
				// if(boundary(mesh, vj.vertex()))
				if(IsBoundaryWithAngleCheck(mesh, vj.vertex()))
				{
					continue;
				}

				HMesh::HalfEdgeID h = vj.opp().halfedge();
				HMesh::VertexID n = vj.vertex();
				double dist = length(mesh, h);

				// collapse edge if allowed and needle is present
				if(dist < low && precond_collapse_edge(mesh, h))
				{
					bool CheckIllegal = true;
					bool illegal = false;

					if (CheckIllegal)
					{
						// Check if the collapse will create a long edge
						HMesh::Walker w = mesh.walker(h); 
						CGLA::Vec3d mid_point = 0.5*(mesh.pos(w.vertex())+mesh.pos(w.opp().vertex()));


						for(HMesh::Walker wc=mesh.walker(w.vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
						{
							if(length(mesh.pos(wc.vertex())-mid_point) > high)
								illegal = true;
						}
						for(HMesh::Walker wc=mesh.walker(w.opp().vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
						{
							if(length(mesh.pos(wc.vertex())-mid_point) > high)
								illegal = true;
						}
					}

					if (!illegal)
					{
						if(vertex_error(mesh, *v, CGLA::Vec3d(mesh.pos(n))) < vertex_error(mesh, n, CGLA::Vec3d(mesh.pos(*v))))
						{
							mesh.pos(*v) = mesh.pos(n);
						}
						mesh.collapse_edge(h);
						did_work = true;
						collapses++;
						break;
					}
					else
					{
						PrevCol++;
					}
				}

			}
		}
	} while(did_work);
//	std::cout << "Collapses: " << collapses << " prevented collapses " << PrevCol << std::endl;
	std::cout << "(" << collapses << "," << PrevCol << ")";
}


void CGELRemeshing::CollapseShortEdgesWithCurvature( HMesh::Manifold &mesh, vtkPointLocator * locator, vtkDoubleArray *scalars )
{
	std::cout << "Collapse";
	// We do not do it exactly as stated as above. We do not check if we create long edges
	//	remove_needles(mesh, low);

	bool did_work = false;
	int collapses = 0;
	int PrevCol = 0;
//	int nRuns = 0;

	// remove needles until no more can be removed
	do
	{
//		int debugc = 0;
//		std::cout << "C";
		did_work = false;
		for(HMesh::VertexIDIterator v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
		{
//			std::cout << " " << debugc++;

			// don't attempt to remove needles if vertex is boundary
			if(boundary(mesh, *v))
			{
				continue;
			}

			for(HMesh::Walker vj = mesh.walker(*v); !vj.full_circle(); vj = vj.circulate_vertex_cw())
			{
				// don't attempt to remove needles if vertex of jumper halfedge is boundary
				if(boundary(mesh, vj.vertex()))
				{
					continue;
				}

				HMesh::HalfEdgeID h = vj.opp().halfedge();
				HMesh::VertexID n = vj.vertex();
				double dist = length(mesh, h);

				CGLA::Vec3d mid_point = 0.5*(mesh.pos(vj.vertex())+mesh.pos(vj.opp().vertex()));
				double p[3];
				p[0] = mid_point[0];
				p[1] = mid_point[1];
				p[2] = mid_point[2];

				int id = locator->FindClosestPoint(p);
				double trgVal = scalars->GetValue(id);
				if (trgVal == 0)
				{
					std::cerr << "Target value of 0 found on target mesh - can not remesh to 0 length. Setting it to 1!!!!!!" << std::endl;
					trgVal = 1;
				}

				double low = 4.0/5.0 * trgVal;
				
				// collapse edge if allowed and needle is present
				if(dist < low && precond_collapse_edge(mesh, h))
				{
					bool CheckIllegal = true;
					bool illegal = false;

					if (CheckIllegal)
					{
						// Check if the collapse will create a long edge
						HMesh::Walker w = mesh.walker(h); 
						CGLA::Vec3d mid_point = 0.5*(mesh.pos(w.vertex())+mesh.pos(w.opp().vertex()));

						//p[0] = mid_point[0];
						//p[1] = mid_point[1];
						//p[2] = mid_point[2];

						//id = locator->FindClosestPoint(p);
						//trgVal = scalars->GetValue(id);
						double high	= 4.0/3.0 * trgVal;  // trgVal value same as before

						for(HMesh::Walker wc=mesh.walker(w.vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
						{
							if(length(mesh.pos(wc.vertex())-mid_point) > high)
								illegal = true;
						}
						for(HMesh::Walker wc=mesh.walker(w.opp().vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
						{
							if(length(mesh.pos(wc.vertex())-mid_point) > high)
								illegal = true;
						}
					}

					if (!illegal)
					{
						if(vertex_error(mesh, *v, CGLA::Vec3d(mesh.pos(n))) < vertex_error(mesh, n, CGLA::Vec3d(mesh.pos(*v))))
						{
							mesh.pos(*v) = mesh.pos(n);
						}
						mesh.collapse_edge(h);
						did_work = true;
						collapses++;
						break;
					}
					else
					{
						PrevCol++;
					}
				}

			}
		}
	} while(did_work);
	//	std::cout << "Collapses: " << collapses << " prevented collapses " << PrevCol << std::endl;
	std::cout << "(" << collapses << "," << PrevCol << ")";
}



void CGELRemeshing::BoundarySmoothing( HMesh::Manifold &mesh )
{
	std::cout << "BoundSmooth";

	std::vector<CGLA::Vec3d> pos(mesh.allocated_vertices());
	std::vector<bool> changed(mesh.allocated_vertices(), false);

	// Run through vertices and check for boundary
	int i = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi, i++)
	{
		//if(boundary(mesh, *vi))
		if(IsBoundaryWithAngleCheck(mesh, *vi))
		{
			std::vector<CGLA::Vec3d> neighbours(2);

			// Test number of neighbours that are also on the boundary
			int n=0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				//if(boundary(mesh, w.vertex()))
				if(IsBoundaryWithAngleCheck(mesh, w.vertex()))
				{
					if (n < 2)
					{
						neighbours[n] = mesh.pos(w.vertex());
					}
					n++;
				}
			}

			// Compute angle between neighbouring edges
			if (n == 2)
			{
				CGLA::Vec3d v0 = mesh.pos(*vi);
				CGLA::Vec3d v1 = neighbours[0];
				CGLA::Vec3d v2 = neighbours[1];
				double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
				double Angle = vtkMath::DegreesFromRadians(a0);

				if (Angle > 170)
				{
					CGLA::Vec3d v = mesh.pos(*vi);

					// Set to average of neighbours
					pos[i] = (v1+v2) / 2;
					changed[i] = true;
				}
			}
		}

	}

	int moved = 0;
	i = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi, i++)
	{
		if (changed[i])
		{
			mesh.pos(*vi) = pos[i];
			moved++;
		}
	}
	std::cout << "(" << moved << ")";
}



void CGELRemeshing::CollapseShortBoundaryEdges( HMesh::Manifold &mesh, double low, double high )
{
	std::cout << "ColBound";

	int collaps = 0;
	int prevents = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
//		if(boundary(mesh, *vi))
		if(IsBoundaryWithAngleCheck(mesh, *vi))
		{
			std::vector<CGLA::Vec3d> neighbours(2);

			int n=0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				// if(boundary(mesh, w.vertex()))
				if(IsBoundaryWithAngleCheck(mesh, w.vertex()))
				{
					if (n < 2)
					{
						neighbours[n] = mesh.pos(w.vertex());
					}
					n++;
				}
			}

			// Compute angle between neighbouring edges
			if (n == 2)
			{
				CGLA::Vec3d v0 = mesh.pos(*vi);
				CGLA::Vec3d v1 = neighbours[0];
				CGLA::Vec3d v2 = neighbours[1];
				double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
				double Angle = vtkMath::DegreesFromRadians(a0);

				// What will the new dist be
				double newdist = length(v1-v2);
				if (newdist > high)
				{
					prevents++;
				}
				else if (Angle > 170)
				{
					HMesh::HalfEdgeID h = w.halfedge();
					double dist = length(mesh, h);
					{
						// collapse edge if allowed
						if(dist < low && precond_collapse_edge(mesh, h))
						{
							mesh.collapse_edge(h);
							collaps++;
						}
					}
				}
			}
		}
	}
	std::cout << "(" << collaps << "," << prevents <<  ")";
}

void CGELRemeshing::CollapseShortBoundaryEdgesWithCurvature(HMesh::Manifold &mesh, vtkPointLocator * locator, vtkDoubleArray *scalars)
{
	std::cout << "ColBound";

	int collaps = 0;
	int prevents = 0;
	for (HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		if (boundary(mesh, *vi))
		{
			std::vector<CGLA::Vec3d> neighbours(2);

			int n = 0;
			HMesh::Walker w = mesh.walker(*vi);
			for (; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				if (boundary(mesh, w.vertex()))
				{
					if (n < 2)
					{
						neighbours[n] = mesh.pos(w.vertex());
					}
					n++;
				}
			}

			// Compute angle between neighbouring edges
			if (n == 2)
			{
				CGLA::Vec3d v0 = mesh.pos(*vi);
				CGLA::Vec3d v1 = neighbours[0];
				CGLA::Vec3d v2 = neighbours[1];
				double a0 = acos(dot(v1 - v0, v2 - v0) / (length(v1 - v0)*length(v2 - v0)));
				double Angle = vtkMath::DegreesFromRadians(a0);
	
				double p[3];
				p[0] = v0[0];
				p[1] = v0[1];
				p[2] = v0[2];

				int id = locator->FindClosestPoint(p);
				double trgVal = scalars->GetValue(id);
				if (trgVal == 0)
				{
					std::cerr << "Target value of 0 found on target mesh - can not remesh to 0 length. Setting it to 1!!!!!!" << std::endl;
					trgVal = 1;
				}


				double low = 4.0 / 5.0 * trgVal;
				double high = 4.0 / 3.0 * trgVal; 

				// What will the new dist be
				double newdist = length(v1 - v2);
				if (newdist > high)
				{
					prevents++;
				}
				else if (Angle > 170)
				{
					HMesh::HalfEdgeID h = w.halfedge();
					double dist = length(mesh, h);
					{
						// collapse edge if allowed
						if (dist < low && precond_collapse_edge(mesh, h))
						{
							mesh.collapse_edge(h);
							collaps++;
						}
						else
						{
							prevents++;
						}
					}
				}
			}
		}
	}
	std::cout << "(" << collaps << "," << prevents << ")";
}


// Test function
//void CGELRemeshing::RemoveCaps(HMesh::Manifold &mesh)
//{
//	std::cout << "removecaps.";
//	double ThresAngle = vtkMath::RadiansFromDegrees(170.0);
//
//	remove_caps(mesh, ThresAngle);
//}
//

/*
The equalize valences() function equalizes the vertex valences by flipping edges. The target
valence target val(v) is 6 and 4 for interior and boundary vertices, respectively. The algorithm
tentatively flips each edge e and checks whether the deviation to the target valences decreases.
If not, the edge is flipped back.

equalize valences()
for each edge e do
let a,b,c,d be the vertices of the two triangles adjacent to e
deviation pre = abs(valence(a)-target val(a)) + abs(valence(b)-target val(b))
+ abs(valence(c)-target val(c)) + abs(valence(d)-target val(d))

flip(e)
deviation post = abs(valence(a)-target val(a)) + abs(valence(b)-target val(b))
+ abs(valence(c)-target val(c)) + abs(valence(d)-target val(d))
if deviation pre <= deviation post do
*/
void CGELRemeshing::EqualizeValences(HMesh::Manifold &mesh )
{
	std::cout << "EqVal";
	//int flips = HMesh::optimize_valency(mesh);
	//std::cout << "(" << flips << ")";
	HMesh::optimize_valency(mesh);
}


//void CGELRemeshing::LaplacianSmooth(HMesh::Manifold& m, double weight, int max_iter)
//{
//	auto vertex_ids = batch_vertices(m);
//	auto new_pos = m.positions_attribute_vector();
//	auto f = [&](const vector<HMesh::VertexID>& vids) {
//		for (HMesh::VertexID v : vids)
//			if (!boundary(m, v))
//				new_pos[v] = m.pos(v) + weight*laplacian(m, v);
//			else
//				new_pos[v] = m.pos(v)
//	};
//
//	for (auto _ : range(0, max_iter)) {
//		for_each_vertex_parallel(CORES, vertex_ids, f);
//		swap(m.positions_attribute_vector(), new_pos);
//	}
//}


void  CGELRemeshing::TAL_smoothing_RAPA_Version(HMesh::Manifold& m, float w, int max_iter)
{
	for (int iter = 0; iter < max_iter; ++iter) {
		HMesh::VertexAttributeVector<float> vertex_areas;
		HMesh::VertexAttributeVector<CGLA::Vec3d> laplacians;

		for (HMesh::VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
		{
			vertex_areas[*vid] = 0;
			for (HMesh::Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
				if (w.face() != HMesh::InvalidFaceID)
					vertex_areas[*vid] += area(m, w.face());
		}

		for (HMesh::VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
		{
			laplacians[*vid] = CGLA::Vec3d(0);
			double weight_sum = 0.0;
			if (IsBoundaryWithAngleCheck(m, *vid))
			{
				double angle_sum = 0;
				for (HMesh::Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
				{
					if (w.face() != HMesh::InvalidFaceID)
					{
						CGLA::Vec3d vec_a = m.pos(w.vertex()) - m.pos(*vid);
						CGLA::Vec3d vec_b = m.pos(w.circulate_vertex_ccw().vertex()) - m.pos(*vid);
						if (vec_a.length() && vec_b.length())
						{
							CGLA::Vec3d vec_a = normalize(m.pos(w.vertex()) - m.pos(*vid));
							CGLA::Vec3d vec_b = normalize(m.pos(w.circulate_vertex_ccw().vertex()) -
								m.pos(*vid));
							angle_sum += acos(std::max(-1.0, std::min(1.0, dot(vec_a, vec_b))));

							//if (std::isnan(angle_sum))
							//{
							//    std::cout << "TAL_smoothing: angle_sum isnan (boundary)" << std::endl;
							//    std::cout << "dot(vec_a,vec_b) " << dot(vec_a, vec_b) << std::endl;
							//    std::cout << "vec a: " << vec_a[0] << vec_a[1] << vec_a[2] << std::endl;
							//    std::cout << "vec b: " << vec_b[0] << vec_b[1] << vec_b[2] << std::endl;
							//}
						}
					}
					if (IsBoundaryWithAngleCheck(m, w.vertex()))
					{
						laplacians[*vid] += m.pos(w.vertex()) - m.pos(*vid);
						weight_sum += 1.0;
					}
				}
				//Vec3d LP1 = laplacians[*vid]; // RASMUS DEBUG

				laplacians[*vid] /= weight_sum;
				//if (!weight_sum)
				//    std::cout << "TAL_smoothing: Weight_sum = 0 (boundary)" << std::endl;
				laplacians[*vid] *= exp(-3.0 * CGLA::sqr(std::max(0.0, M_PI - angle_sum)));

				//if (std::isnan(laplacians[*vid][0]) || std::isnan(laplacians[*vid][1]) || std::isnan(laplacians[*vid][2]))
				//{
				//    std::cout << "TAL_smoothing: laplacians[*vid] isnan (boundary)" << std::endl;
				//    std::cout << "weight_sum: " << weight_sum << std::endl;
				//    std::cout << "angle_sum: " << angle_sum << std::endl;
				//    std::cout << "exp(-3.0*sqr(max(0.0, M_PI-angle_sum)))" << exp(-3.0 * sqr(max(0.0, M_PI - angle_sum))) << std::endl;
				//    std::cout << "laplacian before: " << LP1[0] << LP1[1] << LP1[2] << std::endl;
				//}
			}
			else
			{
				for (HMesh::Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
				{
					float weight = vertex_areas[w.vertex()];
					CGLA::Vec3d l = m.pos(w.vertex()) - m.pos(*vid);
					laplacians[*vid] += weight * l;
					weight_sum += weight;
				}
				//if (!weight_sum)
				//    std::cout << "TAL_smoothing: Weight_sum = 0 (non boundary)" << std::endl;
				laplacians[*vid] /= weight_sum;

				//if (std::isnan(laplacians[*vid][0]) || std::isnan(laplacians[*vid][1]) || std::isnan(laplacians[*vid][2]))
				//{
				//    std::cout << "TAL_smoothing: laplacians[*vid] isnan (not boundary)" << std::endl;
				//    std::cout << "weight_sum: " << weight_sum << std::endl;
				//}


				//                Vec3d n = normal(m, *vid);
				//                if(sqr_length(n)>0.9)
				//                    laplacians[*vid] -= n * dot(n, laplacians[*vid]);
			}
		}
		for (HMesh::VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
		{
			m.pos(*vid) += w * laplacians[*vid];

			//Vec3d LP = laplacians[*vid];
			//Manifold::Vec p = m.pos(*vid);

			//if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2]))
			//{
			//    std::cout << "TAL_smoothing: p isnan" << std::endl;
			//    std::cout << "w: " << w << std::endl;
			//    std::cout << "laplacian: " << LP[0] << LP[1] << LP[2] << std::endl;
			//}
		}
	}
}


/*
The tangential relaxation() function applies an iterative smoothing filter to the mesh. Here
the vertex movement has to be constrained to the vertex' tangent plane in order to stabilize
the following projection operator. Let p be an arbitrary vertex in the current mesh, let n be
its normal, and let q be the position of the vertex as calculated by a smoothing algorithms
with uniform Laplacian weights (see Chapter 7). The new position p0 of p is then computed by
projecting q onto p's tangent plane

	p0 = q + nnT (p - q) :
Again, this can be easily implemented:

tangential relaxation()
for each vertex v do
	q[v] = the barycenter of v's neighbor vertices
for each vertex v do
	let p[v] and n[v] be the position and normal of v, respectively
	p[v] = q[v] + dot(n[v],(p[v]-q[v]))*n[v]
*/
void CGELRemeshing::TangentialRelaxation(vtkImplicitFunction *func, HMesh::Manifold &mesh)
{
	std::cout << "TangRelax.";
	int SmoothType = 2;

	if (SmoothType == 0)
	{
		HMesh::laplacian_smooth(mesh, 0.5, 2); // Shrinks mesh currently (11/11-2015) due to missing boundary check
		return;
	}
	if (SmoothType == 1)
	{
		HMesh::taubin_smooth(mesh);
		return;
	}
	else if (SmoothType == 2)
	{
		HMesh::TAL_smoothing(mesh, 0.5);
		return;
	}

	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi)
	{
		if(!boundary(mesh, *vi))
		{
			// Calculate barycenter (q)
			CGLA::Vec3d q(0);
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				q += mesh.pos(w.vertex());
			}
			if (!w.no_steps())
				std::cout << "TangentialRelaxation w.no_steps = 0" << std::endl;
			q /= w.no_steps();
						
			CGLA::Vec3d p = mesh.pos(*vi);

			// use the implicit function to get a normal
			double pp[3];
			pp[0] = p[0]; pp[1] = p[1]; pp[2] = p[2];

			double *nn = func->FunctionGradient(pp);
			vtkMath::Normalize(nn);

			CGLA::Vec3d norm(-nn[0], -nn[1], -nn[2]);

			CGLA::Vec3d tv = p - q;
		
			CGLA::Vec3d npos = q + CGLA::dot(norm, tv) * norm;

			if (npos.length() > 10000)
			{
				std::cout << "TangentialRelaxation: npos.length > 10000" << std::endl;
			}
			mesh.pos(*vi) = npos;
		}
		//else // RAPA DEBUG 27/1-2019
		//{
		//	CGLA::Vec3d p = mesh.pos(*vi);

		//	if (p.length() > 10000)
		//	{
		//		std::cout << "TangentialRelaxation: p.length > 10000 (boundary)" << std::endl;
		//	}
		//}
	}
}

void CGELRemeshing::TangentialRelaxationOnlyMesh(HMesh::Manifold &mesh)
{
	std::cout << "TangRelax.";
	int SmoothType = 2;

	if (SmoothType == 0)
	{
		HMesh::laplacian_smooth(mesh, 0.5, 2); // Shrinks mesh currently (11/11-2015) due to missing boundary check
		return;
	}
	else if (SmoothType == 1)
	{
		HMesh::taubin_smooth(mesh);
		return;
	}
	else if (SmoothType == 2)
	{
		HMesh::TAL_smoothing(mesh, 0.5);
		return;
	}

	//for (HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	//{
	//	if (!boundary(mesh, *vi))
	//	{
	//		// Calculate barycenter (q)
	//		CGLA::Vec3d q(0);
	//		HMesh::Walker w = mesh.walker(*vi);
	//		for (; !w.full_circle(); w = w.circulate_vertex_ccw())
	//		{
	//			q += mesh.pos(w.vertex());
	//		}
	//		q /= w.no_steps();

	//		CGLA::Vec3d p = mesh.pos(*vi);

	//		// use the implicit function to get a normal
	//		double pp[3];
	//		pp[0] = p[0]; pp[1] = p[1]; pp[2] = p[2];

	//		double *nn = func->FunctionGradient(pp);
	//		vtkMath::Normalize(nn);

	//		CGLA::Vec3d norm(-nn[0], -nn[1], -nn[2]);

	//		CGLA::Vec3d tv = p - q;

	//		CGLA::Vec3d npos = q + CGLA::dot(norm, tv) * norm;

	//		mesh.pos(*vi) = npos;
	//	}
	//}
}


//bool IsNumber(double x)     
//{       
//	// This looks like it should always be true,         
//	// but it's false if x is a NaN.        
//	return (x == x);     
//}   
//
//bool IsFiniteNumber(double x)     
//{        
//	return (x <= DBL_MAX && x >= -DBL_MAX);     
//}        


double CGELRemeshing::VertexArea(HMesh::Manifold &mesh, HMesh::VertexID vi )
{
	double A = 0;
	HMesh::Walker w = mesh.walker(vi);
	for(; !w.full_circle(); w = w.circulate_vertex_ccw())
	{
		if (w.face() != HMesh::InvalidFaceID) // debug check to avoid getting non referenced faces
			A += HMesh::area(mesh, w.face());
	}
	return A;
}

bool CGELRemeshing::IsBoundaryWithAngleCheck(HMesh::Manifold& mesh, HMesh::VertexID vi)
{

	HMesh::Walker w = mesh.walker(vi);
	for (; !w.full_circle(); w = w.circulate_vertex_ccw())
	{
		HMesh::FaceID f1 = w.face();
		if (f1 == HMesh::InvalidFaceID)
		{
			return true;
		}
		HMesh::FaceID f2 = w.opp().face();
		if (f2 == HMesh::InvalidFaceID)
		{
			return true;
		}
		CGLA::Vec3d n1 = HMesh::normal(mesh, f1);
		CGLA::Vec3d n2 = HMesh::normal(mesh, f2);

		// Check if normals points the same way
		if ((n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) < -0.8)
		{
			return true;
		}
	}

	return false;
}



void CGELRemeshing::MoveToVertexCenterMass( vtkImplicitFunction *func, HMesh::Manifold &mesh )
{
	std::cout << "MoveCMS.";

	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi)
	{
		if(!boundary(mesh, *vi))
		{
			CGLA::Vec3d q(0);

			double totArea = 0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				double At = VertexArea(mesh, w.vertex());
				q += mesh.pos(w.vertex()) * At;
				totArea += At;
			}
			if (totArea > 0)
			{
				q = q / (totArea);
				mesh.pos(*vi) = q;
			}
		}
	}
}

void CGELRemeshing::AreaWeigthedTangentialRelaxationOnSphere( HMesh::Manifold &mesh )
{
	std::cout << "AreaTanRel.";

	//	int i = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi)
	{
		if(!boundary(mesh, *vi))
		{
			// Calculate barycenter (q)
			CGLA::Vec3d q(0);
			double totArea = 0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				double At = VertexArea(mesh, w.vertex());
				q += mesh.pos(w.vertex()) * At;
				totArea += At;
			}
			if (totArea > 0)
			{
				q = q / (totArea);
				mesh.pos(*vi) = q;

				CGLA::Vec3d p = mesh.pos(*vi);

				double nn[3];
				nn[0] = p[0]; nn[1] = p[1]; nn[2] = p[2];

				vtkMath::Normalize(nn);

				CGLA::Vec3d norm(-nn[0], -nn[1], -nn[2]);

				CGLA::Vec3d tv = p - q;

				CGLA::Vec3d npos = q + CGLA::dot(norm, tv) * norm;

				mesh.pos(*vi) = npos;
			}
		}
	}
}

void CGELRemeshing::AreaWeigthedTangentialRelaxation( vtkImplicitFunction *func, HMesh::Manifold &mesh )
{
	std::cout << " AreaTanRel.";

	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi)
	{
		if (!boundary(mesh, *vi))
		{
			CGLA::Vec3d q(0);
			double totArea = 0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				double At = VertexArea(mesh, w.vertex());
				q += (mesh.pos(w.vertex()) * At);
				totArea += At;
			}
			if (totArea > 0)
			{
				q = q / totArea;

				CGLA::Vec3d p = mesh.pos(*vi);

				double pp[3];
				pp[0] = p[0]; pp[1] = p[1]; pp[2] = p[2];

				double *nn = func->FunctionGradient(pp);

				vtkMath::Normalize(nn);

				CGLA::Vec3d norm(-nn[0], -nn[1], -nn[2]);
				CGLA::Vec3d tv = p - q;
				CGLA::Vec3d npos = q + CGLA::dot(norm, tv) * norm;

				mesh.pos(*vi) = npos;
			}
		}
	}
}


void CGELRemeshing::BoundaryRelaxation(vtkImplicitFunction *func, HMesh::Manifold &mesh)
{
	// Store the new positions
	std::vector<CGLA::Vec3d> pos(mesh.allocated_vertices());
	std::vector<bool> changed(mesh.allocated_vertices(), false);

	int i=0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi, i++)
	{
		if(boundary(mesh, *vi))
		{
			CGLA::Vec3d q(0);
			int n=0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				if(boundary(mesh, w.vertex()))
				{
					q += mesh.pos(w.vertex());
					++n;
				}
			}
			if (n > 0)
			{
				q = q / n;

				CGLA::Vec3d p = mesh.pos(*vi);

				double pp[3];
				pp[0] = p[0]; pp[1] = p[1]; pp[2] = p[2];

				double *nn = func->FunctionGradient(pp);
				vtkMath::Normalize(nn);

				CGLA::Vec3d norm(-nn[0], -nn[1], -nn[2]);

				CGLA::Vec3d tv = p - q;

				CGLA::Vec3d npos = q + CGLA::dot(norm, tv) * norm;

				pos[i] = npos;
				changed[i] = true;
			}
		}
	}

	i=0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi, i++)
	{
		if (changed[i])
		{
			mesh.pos(*vi) = pos[i];
		}
	}
}


void CGELRemeshing::ProjectToSurfaceOnMesh(HMesh::Manifold &mesh, vtkCellLocator *cLocator, double SanityJump)
{
	bool debug = false;
	std::cout << "ProjSurf";

	for (HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		CGLA::Vec3d pos = mesh.pos(*vi);

		vtkIdType cell_id;
		int sub_id;
		double dist2, totaldist = 0;
		double P[3];
		P[0] = pos[0]; 		
		P[1] = pos[1];
		P[2] = pos[2];

		double tcp[3];
		// see if it is on a border
		cLocator->FindClosestPoint(P, tcp, cell_id, sub_id, dist2);

		bool ok = dist2 < SanityJump * SanityJump;
		if (ok)
		{
			mesh.pos(*vi) = CGLA::Vec3d(tcp[0], tcp[1], tcp[2]);
		}
	}
}

/*
The project to surface() function maps the vertices back to the surface.
*/
void CGELRemeshing::ProjectToSurface( vtkImplicitFunction *func, HMesh::Manifold &mesh, double SanityJump, int iterations)
{
	std::cout << "ProjSurf";

	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi)
	{
		CGLA::Vec3d pos = mesh.pos(*vi);

		double p[3];
		p[0] = pos[0]; p[1] = pos[1];  p[2] = pos[2]; 
		double cp[3];

		// double accuracy = 0.001;
		//double accuracy = 0.0001;
		// double accuracy = SanityJump / 10000;
		double accuracy = SanityJump / 100;  // TODO yet another constant set by trial and error. To avoid sliding along the zero level
		bool ok = CImplicitFunctionUtils::FindClosestPoint(func, p, cp, accuracy, SanityJump, 1e8, iterations);
		if (ok)
		{
			mesh.pos(*vi) = CGLA::Vec3d(cp[0], cp[1], cp[2]);
		}
	}
}

void CGELRemeshing::ProjectToSurfaceByBisection(vtkImplicitFunction* func, HMesh::Manifold& mesh, double epsilon)
{
	std::cout << "ProjBisection";

	for (HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		CGLA::Vec3d pos = mesh.pos(*vi);

		double p[3];
		p[0] = pos[0]; p[1] = pos[1];  p[2] = pos[2];
		double cp[3];

		bool ok = CImplicitFunctionUtils::FindClosestPointByBisection(func, p, cp, epsilon);
		if (ok)
		{
			mesh.pos(*vi) = CGLA::Vec3d(cp[0], cp[1], cp[2]);
		}
	}
}


void CGELRemeshing::ProjectToSurfaceOnSphere( double radius, HMesh::Manifold &mesh, double SanityJump )
{
	std::cout << "projsurf.";

	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi)
	{
		CGLA::Vec3d pos = mesh.pos(*vi);

		double nn[3];
		nn[0] = pos[0]; nn[1] = pos[1]; nn[2] = pos[2];
		vtkMath::Normalize(nn);

		nn[0] *= radius;
		nn[1] *= radius;
		nn[2] *= radius;

		mesh.pos(*vi) = CGLA::Vec3d(nn[0], nn[1], nn[2]);
	}
}


bool CGELRemeshing::VTK2GEL( vtkPolyData *pd, HMesh::Manifold &mesh )
{
//	std::cout << "You are using VTK2GEL : it is strongly recommended to use VTK2GELFaceByFace" << std::endl;
	if (pd->GetNumberOfPoints() == 0)
	{
		std::cerr << "VTK2GEL:  no points in VTK data" << std::endl;
		return false;
	}
	std::vector<CGLA::Vec3d> verts;
	std::vector<int> faces;
	std::vector<int> indices;

	// Extract verts
	for (int i = 0; i < pd->GetNumberOfPoints(); i++)
	{
		double p[3];
		pd->GetPoint(i, p);
		CGLA::Vec3d v(p[0], p[1], p[2]);

		verts.push_back(v);
	}

	// Extract polys
	vtkCellArray *polys = pd->GetPolys();
	polys->InitTraversal();

	vtkIdType npts=0;
	const vtkIdType *pts=0;

	bool stop = false;
	do 
	{
		int res = polys->GetNextCell(npts, pts);

		if (res)
		{
			faces.push_back(npts);
			for (int i = 0; i < npts; i++)
			{
				indices.push_back(pts[i]);
			}
		}
		else
		{
			stop = true;
		}
	} while (!stop);

	mesh.clear();
	//mesh.build(verts.size(), reinterpret_cast<double*>(&verts[0]), faces.size(), &faces[0], &indices[0]);
	HMesh::build(mesh, verts.size(), reinterpret_cast<double*>(&verts[0]), faces.size(), &faces[0], &indices[0]);


	return true;
}


bool CGELRemeshing::VTK2GELFaceByFace(vtkPolyData * pd, HMesh::Manifold & mesh)
{
	if (pd->GetNumberOfPoints() == 0)
	{
		std::cerr << "VTK2GELFaceByFace:  no points in VTK data" << std::endl;
		return false;
	}
	mesh.clear();

	// Extract polys
	vtkCellArray *polys = pd->GetPolys();
	polys->InitTraversal();

	vtkIdType npts = 0;
	const vtkIdType *pts = 0;

	int nfac = 0;
	bool stop = false;
	do
	{
		int res = polys->GetNextCell(npts, pts);

		
		if (res)
		{
			std::vector<HMesh::Manifold::Vec> face;
			for (int i = 0; i < npts; i++)
			{
				double p[3];
				pd->GetPoint(pts[i],p);
				CGLA::Vec3d v(p[0], p[1], p[2]);
				face.push_back(v);
			}
			mesh.add_face(face);
			nfac++;
		}
		else
		{
			stop = true;
		}
	} while (!stop);

	std::cout << "Added " << nfac << " faces" << std::endl;
	double rad = 0.00000001;
	int nval = stitch_mesh(mesh, rad);
	if (nval > 0)
	{
		std::cout << "Constructed GEL mesh with " << nval << " edges that could not be stitched" << std::endl;
	}

	return true;
}

void CGELRemeshing::GEL2VTK( HMesh::Manifold &mesh, vtkPolyData *pd )
{
	vtkPolyData *pdt = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();

	HMesh::VertexAttributeVector<int> vmap((int)mesh.allocated_vertices());

	int k=0;

	for(HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		CGLA::Vec3d v = mesh.pos(*vi);

		double p[3];
		p[0] = v[0];
		p[1] = v[1];
		p[2] = v[2];

		pts->InsertNextPoint(p);

		vmap[*vi] = k++;
	}

	for(HMesh::FaceIDIterator f = mesh.faces_begin(); f != mesh.faces_end(); ++f)
	{        
		std::vector<int> verts;

		for (HMesh::Walker w = mesh.walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
		{
			int vertex_idx = vmap[w.vertex()];			
			verts.push_back(vertex_idx);
		}
		if (verts.size() == 3)
		{
			polys->InsertNextCell(3);

			polys->InsertCellPoint(verts[0]);
			polys->InsertCellPoint(verts[1]);
			polys->InsertCellPoint(verts[2]);
		}
	}

	pdt->SetPoints(pts);
	pts->Delete();
	pdt->SetPolys(polys);
	polys->Delete();

	pd->DeepCopy(pdt);
	pdt->Delete();
}
//
//bool CGELRemeshing::RemeshJAB( vtkImplicitFunction *func, vtkPolyData *pdin, vtkPolyData *pdout )
//{
//	bool debug = false;
//	std::string debugDir   = "E:/data/IMM/Remeshing/RemeshDebugDebug/";
//
//	std::cout << "Starting JAB style remesher" << std::endl;
//
//	HMesh::Manifold m;
//
//	if (!VTK2GELFaceByFace(pdin, m))
//		return false;
//
//	if (HMesh::valid(m))
//	{
//		std::cout << "GEL claims that mesh is not valid" << std::endl;
//		std::cout << "Still trying to remesh" << std::endl;
//		//		return false;
//	}
//
//	if (debug)
//	{
//		obj_save(debugDir+"RemeshDebug.obj", m);
//	}
//	
//	std::cout << "Computing median edge length" << std::endl;
//	std::vector<double> edge_lengths;
//	int n=0;
//	for(HMesh::HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h,++n)
//	{
//		edge_lengths.push_back(length(m, *h));
//	}
//	std::sort(edge_lengths.begin(), edge_lengths.end());
//	double med_len = edge_lengths[n/2];
//	std::cout << "Median edge length: " << med_len << std::endl;
//
//	std::cout << "Computing average length" << std::endl;
//	double avgL = HMesh::average_edge_length(m);
//	std::cout << "Average edge length: " << avgL << std::endl;
//
//	double TargetLen = avgL;
//	
//	for(int iter=0; iter < 10; ++iter)
//	{
//		std::cout << "Iteration " << iter << " ";
//
//		std::vector<HMesh::HalfEdgeID> short_edges, long_edges;
//
//		for(HMesh::HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h)
//		{
//			if(*h < m.walker(*h).opp().halfedge())
//			{
//				double l = HMesh::length(m, *h);
//				if(l > (4.0/3.0) * TargetLen)
//					long_edges.push_back(*h);
//				else if(l< (4.0/5.0) * TargetLen)
//					short_edges.push_back(*h);
//			}
//		}
//		std::cout << "Found " << long_edges.size() << " long edges and " << short_edges.size() << " short edges " << std::endl; 
//
//		std::cout << "Splitting long edges" << std::endl;
//		int splits = 0;
//		for(unsigned int i = 0; i < long_edges.size(); ++i)
//		{
//			if(m.in_use(long_edges[i]))
//			{
//				m.split_edge(long_edges[i]);
//				splits++;
//			}
//		}
//		std::cout << "Splitted " << splits << " long edges" << std::endl;
//
//		if (debug)
////		if (debug && iter == 0)
//		{
//			std::ostringstream ss;
//			ss << debugDir + "RemeshDebug_Step1_it" << iter << "_SplitLongEdges.obj";
//			obj_save(ss.str(), m);
//		}
//
//
//		std::cout << "shortest_edge_triangulate" << std::endl;
//		shortest_edge_triangulate(m);
//
//		if (debug)
//		//		if (debug && iter == 0)
//		{
//			std::ostringstream ss;
//			ss << debugDir + "RemeshDebug_Step2_it" << iter << "_ShortestEdgeTriangulate.obj";
//			obj_save(ss.str(), m);
//		}
//
//
//		int collaps  = 0;
//		for(unsigned int i = 0; i < short_edges.size(); ++i)
//		{
//			if(m.in_use(short_edges[i]) && HMesh::precond_collapse_edge(m, short_edges[i]))
//			{
//				HMesh::Walker w = m.walker(short_edges[i]); 
//				CGLA::Vec3d mid_point = 0.5*(m.pos(w.vertex())+m.pos(w.opp().vertex()));
//
//				bool illegal = false;
//
//				for(HMesh::Walker wc=m.walker(w.vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
//				{
//					if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* TargetLen)
//						illegal = true;
//				}
//				for(HMesh::Walker wc=m.walker(w.opp().vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
//				{
//					if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* TargetLen)
//						illegal = true;
//				}
//				if(!illegal)
//				{
//					m.collapse_edge(short_edges[i],true);
//					collaps++;
//				}
//			}
//		}
//		std::cout << "Collapsed " << collaps << " short edges" << std::endl;
//
//		if (debug)
//		//		if (debug && iter == 0)
//		{
//			std::ostringstream ss;
//			ss << debugDir + "RemeshDebug_Step3_it" << iter << "_CollapsedShortEdges.obj";
//			obj_save(ss.str(), m);
//		}
//		
//		std::cout << "maximize min angle" << std::endl;
//		HMesh::maximize_min_angle(m, 0.98);
//
//		if (debug)
//		//		if (debug && iter == 0)
//		{
//			std::ostringstream ss;
//			ss << debugDir + "RemeshDebug_Step4_it" << iter << "_MaxMinAngle.obj";
//			obj_save(ss.str(), m);
//		}
//
//
//		AreaWeigthedTangentialRelaxation(func, m);
//
//		if (debug)
//		//		if (debug && iter == 0)
//		{
//			std::ostringstream ss;
//			ss << debugDir + "RemeshDebug_Step5_it" << iter << "_AreaWeigthedRelax.obj";
//			obj_save(ss.str(), m);
//		}
//
//		ProjectToSurface(func, m, TargetLen);
//
//		if (debug)
//		//		if (debug && iter == 0)
//		{
//			std::ostringstream ss;
//			ss << debugDir + "RemeshDebug_Step6_it" << iter << "_ProjectToSurface.obj";
//			obj_save(ss.str(), m);
//		}
//
//		m.cleanup();
//	}
//
//	m.cleanup();
//	GEL2VTK(m, pdout);
//
//	return true;
//}


bool CGELRemeshing::RemeshDirectMesh(vtkPolyData *pdin, vtkPolyData *pdout, double targetLength)
{
	CMultiTimer::MultiTimer().Start("Remesh", 1);

	bool debug = false;
	std::string debugDir = "d:/data/IMM/Remeshing/debugout/";

	//RemeshJAB(func, pdin, pdout);
	//return true;

	HMesh::Manifold mesh;

	CMultiTimer::MultiTimer().Start("VTK2GELFaceByFace", 1);
	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	CMultiTimer::MultiTimer().End("VTK2GELFaceByFace");

	if (debug)
	{
		obj_save(debugDir + "RemeshDebug.obj", mesh);
	}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	CMultiTimer::MultiTimer().Start("Computing average length", 1);
	std::cout << "Computing average length" << std::endl;
	double avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;
	CMultiTimer::MultiTimer().End("Computing average length");

	//	double targetLength = avgL * 0.8;
	if (targetLength == 0)  // Only change if no targetlength specified
		targetLength = avgL;
	double loTarget = 4.0 / 5.0   * targetLength;
	double hiTarget = 4.0 / 3.0 * targetLength;

	vtkCellLocator *cLocator = vtkCellLocator::New();
	cLocator->SetDataSet(pdin);
	cLocator->SetNumberOfCellsPerBucket(1);
	cLocator->BuildLocator();

	int saveID = 1000;
	int iterations = 10;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "It#" << i << " ";

		CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
		mesh.cleanup();
		CMultiTimer::MultiTimer().End("mesh.cleanup");

		CMultiTimer::MultiTimer().Start("SplitLongEdgesReengineered", 1);
		SplitLongEdgesReengineered(mesh, hiTarget);
		CMultiTimer::MultiTimer().End("SplitLongEdgesReengineered");

		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step1_it" << i << "_SplitLongEdges.obj";
			obj_save(ss.str(), mesh);
		}

		//shortest_edge_triangulate(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step2_it" << i << "_shortest_edge_triangulate.obj";
		//	obj_save(ss.str(), mesh);
		//}

		CMultiTimer::MultiTimer().Start("CollapseShortEdges", 1);
		CollapseShortEdges(mesh, loTarget, hiTarget);
		CMultiTimer::MultiTimer().End("CollapseShortEdges");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("CollapseShortBoundaryEdges", 1);
		CollapseShortBoundaryEdges(mesh, loTarget, hiTarget);
		CMultiTimer::MultiTimer().End("CollapseShortBoundaryEdges");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3a_it" << i << "_CollapseShortBoundaryEdges.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("EqualizeValences", 1);
		EqualizeValences(mesh);
		CMultiTimer::MultiTimer().End("EqualizeValences");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step4_it" << i << "_EqualizeValences.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("TangentialRelaxation", 1);
		TangentialRelaxationOnlyMesh(mesh);
		CMultiTimer::MultiTimer().End("TangentialRelaxation");
		//		AreaWeigthedTangentialRelaxation(func, mesh);
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxation.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("BoundarySmoothing", 1);
		BoundarySmoothing(mesh);
		CMultiTimer::MultiTimer().End("BoundarySmoothing");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5a_it" << i << "_BoundarySmoothing.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("MaximiseMinAngle", 1);
		MaximiseMinAngle(mesh);
		CMultiTimer::MultiTimer().End("MaximiseMinAngle");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step6_it" << i << "_MaximiseMinAngle.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("ProjectToSurface", 1);
		double maxJump = 2 * avgL;
		ProjectToSurfaceOnMesh(mesh, cLocator, maxJump);
		CMultiTimer::MultiTimer().End("ProjectToSurface");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step7_it" << i << "_ProjectToSurface.obj";
			obj_save(ss.str(), mesh);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
	mesh.cleanup();
	CMultiTimer::MultiTimer().End("mesh.cleanup");

	CMultiTimer::MultiTimer().Start("Computing average length", 1);
	avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;
	CMultiTimer::MultiTimer().End("Computing average length");

	//CMultiTimer::MultiTimer().Start("VisualiseBadFaceNormals", 1);
	//int badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
	//std::cout << "Number of bad normals " << badNormals << std::endl;
	//CMultiTimer::MultiTimer().End("VisualiseBadFaceNormals");

	CMultiTimer::MultiTimer().Start("GEL2VTK", 1);
	GEL2VTK(mesh, pdout);
	CMultiTimer::MultiTimer().End("GEL2VTK");

	CMultiTimer::MultiTimer().End("Remesh");
	return true;


}


/*
remesh(target edge length)
low = 4/5 * target edge length
high = 4/3 * target edge length
for i = 0 to 10 do
split long edges(high)
collapse short edges(low,high)
equalize valences()
tangential relaxation()
project to surface()
*/

//#include <HMesh/obj_save.h>

bool CGELRemeshing::Remesh( vtkImplicitFunction *func, vtkPolyData *pdin, vtkPolyData *pdout, double targetEdgeLength, double edgeFactor)
{
	CMultiTimer::MultiTimer().Start("Remesh", 1);

	bool debug = false;
	std::string debugDir   = "d:/data/test/debugout/";

	//RemeshJAB(func, pdin, pdout);
	//return true;

	HMesh::Manifold mesh;

	CMultiTimer::MultiTimer().Start("VTK2GELFaceByFace", 1);
	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	CMultiTimer::MultiTimer().End("VTK2GELFaceByFace");

	//CMultiTimer::MultiTimer().Start("VTK2GEL", 1);
	//std::cout << "VTK2GEL" << std::endl;
	//if (!VTK2GEL(pdin, mesh))
	//	return false;

	CMultiTimer::MultiTimer().End("VTK2GEL");

	if (debug)
	{
		obj_save(debugDir+"RemeshDebug.obj", mesh);
	}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
//		return false;
	}

	double avgL = targetEdgeLength;
	if (targetEdgeLength <= 0)
	{
		CMultiTimer::MultiTimer().Start("Computing average length", 1);
		std::cout << "Computing average length" << std::endl;
		avgL = HMesh::average_edge_length(mesh);
		std::cout << "Average edge length: " << avgL << std::endl;
		CMultiTimer::MultiTimer().End("Computing average length");
	}

//	double targetLength = avgL * 0.8;
	double targetLength = avgL * edgeFactor;
	double loTarget = 4.0 / 5.0 * targetLength;
	double hiTarget = 4.0 / 3.0 * targetLength;

	int saveID = 1000;
	int iterations = 10;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "It#" << i << " ";

		CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
		mesh.cleanup();
		CMultiTimer::MultiTimer().End("mesh.cleanup");

		CMultiTimer::MultiTimer().Start("SplitLongEdgesReengineered", 1);
		SplitLongEdgesReengineered(mesh, hiTarget);
		CMultiTimer::MultiTimer().End("SplitLongEdgesReengineered");

		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step1_it" << i << "_SplitLongEdges.obj";
			obj_save(ss.str(), mesh);
		}

		//shortest_edge_triangulate(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step2_it" << i << "_shortest_edge_triangulate.obj";
		//	obj_save(ss.str(), mesh);
		//}

		CMultiTimer::MultiTimer().Start("CollapseShortEdges", 1);
		CollapseShortEdges(mesh, loTarget, hiTarget);
		CMultiTimer::MultiTimer().End("CollapseShortEdges");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("CollapseShortBoundaryEdges", 1);
		CollapseShortBoundaryEdges(mesh, loTarget, hiTarget);
		CMultiTimer::MultiTimer().End("CollapseShortBoundaryEdges");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3a_it" << i << "_CollapseShortBoundaryEdges.obj";
			obj_save(ss.str(), mesh);
		}
		
		CMultiTimer::MultiTimer().Start("EqualizeValences", 1);
		EqualizeValences(mesh);
		CMultiTimer::MultiTimer().End("EqualizeValences");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step4_it" << i << "_EqualizeValences.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("TangentialRelaxation", 1);
		TangentialRelaxation(func, mesh);
		CMultiTimer::MultiTimer().End("TangentialRelaxation");
//		AreaWeigthedTangentialRelaxation(func, mesh);
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxation.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("BoundarySmoothing", 1);
		BoundarySmoothing(mesh);
		CMultiTimer::MultiTimer().End("BoundarySmoothing");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5a_it" << i << "_BoundarySmoothing.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("MaximiseMinAngle", 1);
		MaximiseMinAngle(mesh);
		CMultiTimer::MultiTimer().End("MaximiseMinAngle");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step6_it" << i << "_MaximiseMinAngle.obj";
			obj_save(ss.str(), mesh);
		}

		int projectionIts = 100;
		double sanityJump = targetLength * 10;
		CMultiTimer::MultiTimer().Start("ProjectToSurface", 1);
		ProjectToSurface(func, mesh, sanityJump, projectionIts);
		CMultiTimer::MultiTimer().End("ProjectToSurface");

		//CMultiTimer::MultiTimer().Start("ProjectToSurface", 1);
		//ProjectToSurface(func, mesh, targetLength);
		//CMultiTimer::MultiTimer().End("ProjectToSurface");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step7_it" << i << "_ProjectToSurface.obj";
			obj_save(ss.str(), mesh);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
	mesh.cleanup();
	CMultiTimer::MultiTimer().End("mesh.cleanup");

	CMultiTimer::MultiTimer().Start("Computing average length", 1);
	avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;
	CMultiTimer::MultiTimer().End("Computing average length");

	//CMultiTimer::MultiTimer().Start("VisualiseBadFaceNormals", 1);
	//int badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
	//std::cout << "Number of bad normals " << badNormals << std::endl;
	//CMultiTimer::MultiTimer().End("VisualiseBadFaceNormals");

	CMultiTimer::MultiTimer().Start("GEL2VTK", 1);
	GEL2VTK(mesh, pdout);
	CMultiTimer::MultiTimer().End("GEL2VTK");

	CMultiTimer::MultiTimer().End("Remesh");
	return true;
}



bool CGELRemeshing::ClosestPointStatistics(vtkPolyData *pd, double &avgDist)
{
	std::cout << "Computing point cloud statistics" << std::endl;

	const int NPoint = pd->GetNumberOfPoints();
	if (!NPoint)
		return false;

	std::vector<double> dists;

	vtkPointLocator* locator = vtkPointLocator::New();
	locator->SetDataSet(pd);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	int Ndists = 0;
	double sumDists = 0;
	for (int i = 0; i < NPoint; i++)
	{
		double p[3];
		pd->GetPoint(i, p);

		vtkIdList* neighPts = vtkIdList::New();
		locator->FindClosestNPoints(2, p, neighPts);

		if (neighPts->GetNumberOfIds() == 2)
		{
			vtkIdType cid = neighPts->GetId(1);
			double cp[3];

			pd->GetPoint(cid, cp);

			double dist = sqrt(vtkMath::Distance2BetweenPoints(p, cp));
			if (dist > 0)
			{
				dists.push_back(dist);
				Ndists++;
				sumDists += dist;
			}
		}
		neighPts->Delete();
	}
	locator->Delete();

	avgDist = 0;
	if (Ndists)
		avgDist = sumDists / double(Ndists);

	return true;
}


bool CGELRemeshing::RemeshWithMultiPointProjection(vtkImplicitFunction* func, vtkPolyData* pdin, vtkPolyData* pdout, double targetEdgeLength)
{
	CMultiTimer::MultiTimer().Start("Remesh", 1);

	bool debug = false;
	std::string debugDir = "d:/data/test/debugout/";

	int NPoints = pdin->GetNumberOfPoints();

	double avgDist = 0;
	ClosestPointStatistics(pdin, avgDist);
	double tolerance = avgDist / 2;
		
	int NTargetPoints = 100000;
	// Start by creating a lot of points

	int NPointPerVertex = NTargetPoints / NPoints;

	vtkPolyData *pd = vtkPolyData::New();
	vtkPoints *newPoints = vtkPoints::New();
	vtkCellArray* verts = vtkCellArray::New();

	double RandSpread = avgDist * 2;
	for (int i = 0; i < pdin->GetNumberOfPoints(); i++)
	{
		double P[3];
		pdin->GetPoint(i, P);

		for (int j = 0; j < NPointPerVertex; j++)
		{
			double Pnew[3];
			Pnew[0] = P[0] + vtkMath::Random(-RandSpread, RandSpread);
			Pnew[1] = P[1] + vtkMath::Random(-RandSpread, RandSpread);
			Pnew[2] = P[2] + vtkMath::Random(-RandSpread, RandSpread);

			// Project point
			double dist = func->FunctionValue(Pnew);
			// First check if we are within reasonable distance
			if (std::abs(dist) < 4 * avgDist)
			{
				double Pn2[3];
				Pn2[0] = Pnew[0];
				Pn2[1] = Pnew[1];
				Pn2[2] = Pnew[2];

				double minDist = dist;
				double minP[3];
				minP[0] = Pnew[0];
				minP[1] = Pnew[1];
				minP[2] = Pnew[2];

				for (int it = 0; it < 10; it++)
				{
					dist = func->FunctionValue(Pn2);
					double* n = func->FunctionGradient(Pn2);
					double nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
					if (nl)
					{
						n[0] /= nl;
						n[1] /= nl;
						n[2] /= nl;

						Pn2[0] = Pn2[0] + dist * n[0];
						Pn2[1] = Pn2[1] + dist * n[1];
						Pn2[2] = Pn2[2] + dist * n[2];
						double newD = func->FunctionValue(Pn2);
						if (newD < minDist)
						{
							minDist = newD;
							minP[0] = Pn2[0];
							minP[1] = Pn2[1];
							minP[2] = Pn2[2];
						}
					}
				}

				double newDist = func->FunctionValue(minP);
				if (newDist < tolerance)
				{
					vtkIdType id = newPoints->InsertNextPoint(minP);
					verts->InsertNextCell(1);
					verts->InsertCellPoint(id);
				}
			}
		}
	}


	//// Project points to zero level
	//for (int i = 0; i < newPoints->GetNumberOfPoints(); i++)
	//{
	//	double P[3];
	//	newPoints->GetPoint(i, P);

	//	double *n = func->FunctionGradient(P);
	//	double nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	//	if (nl)
	//	{
	//		n[0] /= nl;
	//		n[1] /= nl;
	//		n[2] /= nl;

	//		double dist = func->FunctionValue(P);
	//		if (std::abs(dist) < 1e100)
	//		{
	//			double Pnew[3];
	//			Pnew[0] = P[0] + dist * n[0];
	//			Pnew[1] = P[1] + dist * n[1];
	//			Pnew[2] = P[2] + dist * n[2];

	//			newPoints->SetPoint(i, Pnew);
	//		}
	//	}
	//}

	std::cout << "Created " << newPoints->GetNumberOfPoints() << " points from " << pdin->GetNumberOfPoints() << " input points" << std::endl;

	pd->SetPoints(newPoints);
	pd->SetVerts(verts);
	newPoints->Delete();
	verts->Delete();


	if (debug)
	{
		vtkExtMisc::WritePDVTK(pd, debugDir + "MultipointSamples.vtk");
	}

	pdout->DeepCopy(pd);
	pd->Delete();


	//HMesh::Manifold mesh;

	//CMultiTimer::MultiTimer().Start("VTK2GELFaceByFace", 1);
	//std::cout << "VTK2GELFaceByFace" << std::endl;
	//if (!VTK2GELFaceByFace(pdin, mesh))
	//	return false;

	//CMultiTimer::MultiTimer().End("VTK2GELFaceByFace");

	//if (debug)
	//{
	//	obj_save(debugDir + "RemeshDebug.obj", mesh);
	//}

	//if (!valid(mesh))
	//{
	//	std::cout << "GEL claims that mesh is not valid" << std::endl;
	//	std::cout << "Still trying to remesh" << std::endl;
	//	//		return false;
	//}

	//double avgL = targetEdgeLength;
	//if (targetEdgeLength <= 0)
	//{
	//	CMultiTimer::MultiTimer().Start("Computing average length", 1);
	//	std::cout << "Computing average length" << std::endl;
	//	avgL = HMesh::average_edge_length(mesh);
	//	std::cout << "Average edge length: " << avgL << std::endl;
	//	CMultiTimer::MultiTimer().End("Computing average length");
	//}

	////	double targetLength = avgL * 0.8;
	//double targetLength = avgL;
	//double loTarget = 4.0 / 5.0 * targetLength;
	//double hiTarget = 4.0 / 3.0 * targetLength;





	//int saveID = 1000;
	//int iterations = 10;
	//for (int i = 0; i < iterations; i++)
	//{
	//	std::cout << "It#" << i << " ";

	//	CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
	//	mesh.cleanup();
	//	CMultiTimer::MultiTimer().End("mesh.cleanup");

	//	double epsilon = targetLength / 100;
	//	CMultiTimer::MultiTimer().Start("ProjectToSurface", 1);
	//	ProjectToSurfaceByBisection(func, mesh, epsilon);
	//	CMultiTimer::MultiTimer().End("ProjectToSurface");

	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step0_it" << i << "_ProjectToSurface.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_it" << i << "_AngleBoundaries.vtk";
	//		VisualiseAngleBasedBoundary(mesh, ss.str(), true);
	//	}


	//	CMultiTimer::MultiTimer().Start("SplitLongEdgesReengineered", 1);
	//	SplitLongEdgesReengineered(mesh, hiTarget);
	//	CMultiTimer::MultiTimer().End("SplitLongEdgesReengineered");

	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step1_it" << i << "_SplitLongEdges.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	CMultiTimer::MultiTimer().Start("CollapseShortEdges", 1);
	//	CollapseShortEdges(mesh, loTarget, hiTarget);
	//	CMultiTimer::MultiTimer().End("CollapseShortEdges");
	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	CMultiTimer::MultiTimer().Start("CollapseShortBoundaryEdges", 1);
	//	CollapseShortBoundaryEdges(mesh, loTarget, hiTarget);
	//	CMultiTimer::MultiTimer().End("CollapseShortBoundaryEdges");
	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3a_it" << i << "_CollapseShortBoundaryEdges.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	CMultiTimer::MultiTimer().Start("EqualizeValences", 1);
	//	EqualizeValences(mesh);
	//	CMultiTimer::MultiTimer().End("EqualizeValences");
	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step4_it" << i << "_EqualizeValences.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	CMultiTimer::MultiTimer().Start("TangentialRelaxation", 1);
	//	TAL_smoothing_RAPA_Version(mesh, 0.5, 1);
	//	// TangentialRelaxation(func, mesh);
	//	CMultiTimer::MultiTimer().End("TangentialRelaxation");
	//	//		AreaWeigthedTangentialRelaxation(func, mesh);
	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxation.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	CMultiTimer::MultiTimer().Start("BoundarySmoothing", 1);
	//	BoundarySmoothing(mesh);
	//	CMultiTimer::MultiTimer().End("BoundarySmoothing");
	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5a_it" << i << "_BoundarySmoothing.obj";
	//		obj_save(ss.str(), mesh);
	//	}

	//	CMultiTimer::MultiTimer().Start("MaximiseMinAngle", 1);
	//	MaximiseMinAngle(mesh);
	//	CMultiTimer::MultiTimer().End("MaximiseMinAngle");
	//	if (debug)
	//	{
	//		std::ostringstream ss;
	//		ss << debugDir + "RemeshDebug_" << saveID++ << "_Step6_it" << i << "_MaximiseMinAngle.obj";
	//		obj_save(ss.str(), mesh);
	//	}


	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	//CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
	//mesh.cleanup();
	//CMultiTimer::MultiTimer().End("mesh.cleanup");

	//CMultiTimer::MultiTimer().Start("Computing average length", 1);
	//avgL = HMesh::average_edge_length(mesh);
	//std::cout << "Average edge length: " << avgL << std::endl;
	//CMultiTimer::MultiTimer().End("Computing average length");

	//CMultiTimer::MultiTimer().Start("VisualiseBadFaceNormals", 1);
	//int badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
	//std::cout << "Number of bad normals " << badNormals << std::endl;
	//CMultiTimer::MultiTimer().End("VisualiseBadFaceNormals");

	//CMultiTimer::MultiTimer().Start("GEL2VTK", 1);
	//GEL2VTK(mesh, pdout);
	//CMultiTimer::MultiTimer().End("GEL2VTK");

	//CMultiTimer::MultiTimer().End("Remesh");
	return true;
}


bool CGELRemeshing::RemeshWithBisectionProjection(vtkImplicitFunction* func, vtkPolyData* pdin, vtkPolyData* pdout, double targetEdgeLength)
{
	CMultiTimer::MultiTimer().Start("Remesh", 1);

	bool debug = false;
	std::string debugDir = "d:/data/test/debugout/";

	//RemeshJAB(func, pdin, pdout);
	//return true;

	HMesh::Manifold mesh;

	CMultiTimer::MultiTimer().Start("VTK2GELFaceByFace", 1);
	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	CMultiTimer::MultiTimer().End("VTK2GELFaceByFace");

	if (debug)
	{
		obj_save(debugDir + "RemeshDebug.obj", mesh);
	}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	double avgL = targetEdgeLength;
	if (targetEdgeLength <= 0)
	{
		CMultiTimer::MultiTimer().Start("Computing average length", 1);
		std::cout << "Computing average length" << std::endl;
		avgL = HMesh::average_edge_length(mesh);
		std::cout << "Average edge length: " << avgL << std::endl;
		CMultiTimer::MultiTimer().End("Computing average length");
	}

	//	double targetLength = avgL * 0.8;
	double targetLength = avgL;
	double loTarget = 4.0 / 5.0 * targetLength;
	double hiTarget = 4.0 / 3.0 * targetLength;

	int saveID = 1000;
	int iterations = 10;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "It#" << i << " ";

		CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
		mesh.cleanup();
		CMultiTimer::MultiTimer().End("mesh.cleanup");

		double epsilon = targetLength / 100;
		CMultiTimer::MultiTimer().Start("ProjectToSurface", 1);
		ProjectToSurfaceByBisection(func, mesh, epsilon);
		CMultiTimer::MultiTimer().End("ProjectToSurface");

		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step0_it" << i << "_ProjectToSurface.obj";
			obj_save(ss.str(), mesh);
		}

		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_it" << i << "_AngleBoundaries.vtk";
		//	VisualiseAngleBasedBoundary(mesh, ss.str(), true);
		//}


		CMultiTimer::MultiTimer().Start("SplitLongEdgesReengineered", 1);
		SplitLongEdgesReengineered(mesh, hiTarget);
		CMultiTimer::MultiTimer().End("SplitLongEdgesReengineered");

		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step1_it" << i << "_SplitLongEdges.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("CollapseShortEdges", 1);
		CollapseShortEdges(mesh, loTarget, hiTarget);
		CMultiTimer::MultiTimer().End("CollapseShortEdges");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("CollapseShortBoundaryEdges", 1);
		CollapseShortBoundaryEdges(mesh, loTarget, hiTarget);
		CMultiTimer::MultiTimer().End("CollapseShortBoundaryEdges");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3a_it" << i << "_CollapseShortBoundaryEdges.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("EqualizeValences", 1);
		EqualizeValences(mesh);
		CMultiTimer::MultiTimer().End("EqualizeValences");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step4_it" << i << "_EqualizeValences.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("TangentialRelaxation", 1);
		TAL_smoothing_RAPA_Version(mesh, 0.5, 1);
		// TangentialRelaxation(func, mesh);
		CMultiTimer::MultiTimer().End("TangentialRelaxation");
		//		AreaWeigthedTangentialRelaxation(func, mesh);
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxation.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("BoundarySmoothing", 1);
		BoundarySmoothing(mesh);
		CMultiTimer::MultiTimer().End("BoundarySmoothing");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5a_it" << i << "_BoundarySmoothing.obj";
			obj_save(ss.str(), mesh);
		}

		CMultiTimer::MultiTimer().Start("MaximiseMinAngle", 1);
		MaximiseMinAngle(mesh);
		CMultiTimer::MultiTimer().End("MaximiseMinAngle");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step6_it" << i << "_MaximiseMinAngle.obj";
			obj_save(ss.str(), mesh);
		}

		
		std::cout << std::endl;
	}
	std::cout << std::endl;

	CMultiTimer::MultiTimer().Start("mesh.cleanup", 1);
	mesh.cleanup();
	CMultiTimer::MultiTimer().End("mesh.cleanup");

	CMultiTimer::MultiTimer().Start("Computing average length", 1);
	avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;
	CMultiTimer::MultiTimer().End("Computing average length");

	//CMultiTimer::MultiTimer().Start("VisualiseBadFaceNormals", 1);
	//int badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
	//std::cout << "Number of bad normals " << badNormals << std::endl;
	//CMultiTimer::MultiTimer().End("VisualiseBadFaceNormals");

	CMultiTimer::MultiTimer().Start("GEL2VTK", 1);
	GEL2VTK(mesh, pdout);
	CMultiTimer::MultiTimer().End("GEL2VTK");

	CMultiTimer::MultiTimer().End("Remesh");
	return true;
}


bool CGELRemeshing::ProjectOnly(vtkImplicitFunction* func, vtkPolyData* pdin, vtkPolyData* pdout)
{
	CMultiTimer::MultiTimer().Start("ProjectOnly", 1);

	bool debug = false;
	std::string debugDir = "d:/data/test/debugout/";

	HMesh::Manifold mesh;

	CMultiTimer::MultiTimer().Start("VTK2GELFaceByFace", 1);
	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	CMultiTimer::MultiTimer().End("VTK2GELFaceByFace");

	double avgL = 0;
	if (avgL <= 0)
	{
		CMultiTimer::MultiTimer().Start("Computing average length", 1);
		std::cout << "Computing average length" << std::endl;
		avgL = HMesh::average_edge_length(mesh);
		std::cout << "Average edge length: " << avgL << std::endl;
		CMultiTimer::MultiTimer().End("Computing average length");
	}

	double targetLength = avgL;

	if (debug)
	{
		obj_save(debugDir + "RemeshDebug.obj", mesh);
	}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	if (debug)
	{
		std::string oname = debugDir + "MeshGradients.vtk";
		VisualiseGradients(mesh, func, oname);
	}
	

	int saveID = 1000;
	int iterations = 1;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "It#" << i << " ";

		int projectionIts = 100;
		double sanityJump = targetLength * 10;
		double epsilon = targetLength / 100;

		CMultiTimer::MultiTimer().Start("ProjectToSurface", 1);
		// ProjectToSurface(func, mesh, sanityJump, projectionIts);
		ProjectToSurfaceByBisection(func, mesh, epsilon);
		//ProjectToSurfaceByBisection(func, mesh, epsilon);
		//ProjectToSurfaceByBisection(func, mesh, epsilon);
		CMultiTimer::MultiTimer().End("ProjectToSurface");
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step7_it" << i << "_ProjectToSurface.obj";
			obj_save(ss.str(), mesh);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	if (debug)
	{
		std::string oname = debugDir + "MeshGradients_after.vtk";
		VisualiseGradients(mesh, func, oname);
	}


	CMultiTimer::MultiTimer().Start("Computing average length", 1);
	avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;
	CMultiTimer::MultiTimer().End("Computing average length");

	CMultiTimer::MultiTimer().Start("GEL2VTK", 1);
	GEL2VTK(mesh, pdout);
	CMultiTimer::MultiTimer().End("GEL2VTK");

	CMultiTimer::MultiTimer().End("ProjectOnly");
	return true;
}


bool CGELRemeshing::RemeshWithTargetLengthsDirectlyOnMesh(vtkPolyData *pdin, vtkPolyData *curv, vtkPolyData *pdout)
{
	// First check if there are scalars
	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(curv->GetPointData()->GetScalars());
	if (!scalars)
	{
		std::cerr << "No scalars to use for target length" << std::endl;
		return false;
	}

	bool debug = false;
//	bool debug = true;
//	std::string debugDir = "d:/data/IMM/Remeshing/debugout/";
	std::string debugDir = "D:\\data\\IMM\\HRTF Simul\\Fabian13062016\\RemeshDebug\\";

	//RemeshJAB(func, pdin, pdout);
	//return true;

	HMesh::Manifold mesh;

	std::cout << "VTK2GELFaceByFace" << std::endl;
	//if (!VTK2GEL(pdin, mesh))
	//	return false;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;


	//if (debug)
	//{
	//	obj_save(debugDir + "RemeshDebug.obj", mesh);
	//	VisualiseFaceNormals(mesh, debugDir + "RemeshDebug_FaceNormals.vtk");
	//	int badNormals = VisualiseBadFaceNormals(mesh, debugDir + "RemeshDebug_BadFaceNormals.vtk", true);
	//	std::cout << "Bad normals in input mesh " << badNormals << std::endl;
	//}
	//else
	//{
	//	std::cout << "Bad normals in input mesh " << VisualiseBadFaceNormals(mesh, std::string(""), false) << std::endl;
	//}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	// USed for sanity check in projection length jump
	double avgL = HMesh::average_edge_length(mesh);


	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(curv);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	vtkCellLocator *cLocator = vtkCellLocator::New();
	cLocator->SetDataSet(pdin);
	cLocator->SetNumberOfCellsPerBucket(1); 
	cLocator->BuildLocator();

	int saveID = 1000;
	int iterations = 10;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "It#" << i << " ";

		mesh.cleanup();

		int badNormals = 0;

//		int badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);

		//		SplitLongEdgesReengineered(mesh, hiTarget);
		SplitLongEdgesWithCurvature(mesh, locator, scalars);
		if (debug)
		{
			std::ostringstream ss;
			ss << debugDir + "RemeshDebug_" << saveID++ << "_Step1_it" << i << "_SplitLongEdges.obj";
			obj_save(ss.str(), mesh);
//			badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		}

		//shortest_edge_triangulate(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step2_it" << i << "_shortest_edge_triangulate.obj";
		//	obj_save(ss.str(), mesh);
		//}

		//		CollapseShortEdges(mesh, loTarget, hiTarget);
		CollapseShortEdgesWithCurvature(mesh, locator, scalars);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges.obj";
		//	obj_save(ss.str(), mesh);
		//	ss.str("");
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges_badnormals.vtk";
		//	badNormals = VisualiseBadFaceNormals(mesh, ss.str(), true);
		//}

		CollapseShortBoundaryEdgesWithCurvature(mesh, locator, scalars);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3a_it" << i << "_CollapseShortBoundaryEdges.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//	std::cout << "Badnormals after CollapseShortBoundaryEdgesWithCurvature " << badNormals << std::endl;
		//}

		// Dont do this in the first iterations
		if (i > iterations / 2)
		{
			EqualizeValences(mesh);
			//if (debug)
			//{
			//	std::ostringstream ss;
			//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step4_it" << i << "_EqualizeValences.obj";
			//	obj_save(ss.str(), mesh);
			//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
			//}
		}

		TangentialRelaxationOnlyMesh( mesh);
		//		AreaWeigthedTangentialRelaxation(func, mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxation.obj";
		//	obj_save(ss.str(), mesh);
		//	ss.str("");
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxations_badnormals.vtk";
		//	badNormals = VisualiseBadFaceNormals(mesh, ss.str(), true);
		//	std::cout << "Badnormals after TangentialRelaxation " << badNormals << std::endl;
		//}

		BoundarySmoothing(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5a_it" << i << "_BoundarySmoothing.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//	std::cout << "Badnormals after BoundarySmoothing " << badNormals << std::endl;
		//}

		MaximiseMinAngle(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step6_it" << i << "_MaximiseMinAngle.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//	std::cout << "Badnormals after MaximiseMinAngle " << badNormals << std::endl;

		//}

		double maxJump = 2 * avgL;
		ProjectToSurfaceOnMesh(mesh, cLocator, maxJump);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step7_it" << i << "_ProjectToSurface.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//}

		//std::cout << "BadNorms " << VisualiseBadFaceNormals(mesh, std::string(""), false);
		std::cout << std::endl;
	}
	std::cout << std::endl;

	locator->Delete();
	cLocator->Delete();

//	std::cout << "Cleanup" << std::endl;
	mesh.cleanup();
//	std::cout << "Cleanup end" << std::endl;

	//if (debug)
	//{
	//	VisualiseFaceNormals(mesh, debugDir + "RemeshDebug_Remeshed_FaceNormals.vtk");
	//	int badNormals = VisualiseBadFaceNormals(mesh, debugDir + "RemeshDebug_Remeshed_BadFaceNormals.vtk", true);
	//	std::cout << "Number of bad normals " << badNormals << std::endl;
	//}

	//avgL = HMesh::average_edge_length(mesh);
	//std::cout << "Average edge length: " << avgL << std::endl;

//	std::cout << "GEL2VTK" << std::endl;
	GEL2VTK(mesh, pdout);
//	std::cout << "GEL2VTKend" << std::endl;
	return true;

}


bool CGELRemeshing::RemeshWithCurvature( vtkImplicitFunction *func, vtkPolyData *pdin, vtkPolyData *curv, vtkPolyData *pdout )
{
	//bool debug = false;
	//std::string debugDir   = "E:/data/IMM/Remeshing/RemeshDebugDebug/";

	// First check if there are scalars
	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(curv->GetPointData()->GetScalars());
	if (!scalars)
	{
		std::cerr << "No scalars to use for target length" << std::endl;
		return false;
	}

	bool debug = false;
//	bool debug = true;
	std::string debugDir = "d:/data/IMM/Remeshing/debugout/";

	//RemeshJAB(func, pdin, pdout);
	//return true;

	HMesh::Manifold mesh;

	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	//std::cout << "VTK2GEL" << std::endl;
	//if (!VTK2GEL(pdin, mesh))
	//	return false;

	//if (debug)
	//{
	//	obj_save(debugDir+"RemeshDebug.obj", mesh);
	//	VisualiseFaceNormals(mesh, debugDir + "RemeshDebug_FaceNormals.vtk");
	//	int badNormals = VisualiseBadFaceNormals(mesh, debugDir + "RemeshDebug_BadFaceNormals.vtk", true);
	//	std::cout << "Bad normals in input mesh " << badNormals << std::endl;
	//}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	std::cout << "Computing average length" << std::endl;
	double avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;

//	double targetLength = avgL * 0.8;
	double targetLength = avgL;
	double loTarget = 4.0 / 5.0   * targetLength;
	double hiTarget = 4.0 / 3.0 * targetLength;


	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(curv);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();



	int saveID = 1000;
	int iterations = 10;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "It#" << i << " ";

		mesh.cleanup();

		//int badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);

//		SplitLongEdgesReengineered(mesh, hiTarget);
		SplitLongEdgesWithCurvature(mesh, locator, scalars);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step1_it" << i << "_SplitLongEdges.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//}

		//shortest_edge_triangulate(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step2_it" << i << "_shortest_edge_triangulate.obj";
		//	obj_save(ss.str(), mesh);
		//}

//		CollapseShortEdges(mesh, loTarget, hiTarget);
		CollapseShortEdgesWithCurvature(mesh, locator, scalars);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges.obj";
		//	obj_save(ss.str(), mesh);
		//	ss.str("");
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3_it" << i << "_CollapseShortEdges_badnormals.vtk";
		//	badNormals = VisualiseBadFaceNormals(mesh, ss.str(), true);
		//}

		CollapseShortBoundaryEdgesWithCurvature(mesh, locator, scalars);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step3a_it" << i << "_CollapseShortBoundaryEdges.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//	std::cout << "Badnormals after CollapseShortBoundaryEdgesWithCurvature " << badNormals << std::endl;
		//}

		// Dont do this in the first iterations
		if (i > iterations / 2)
		{
			EqualizeValences(mesh);
			//if (debug)
			//{
			//	std::ostringstream ss;
			//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step4_it" << i << "_EqualizeValences.obj";
			//	obj_save(ss.str(), mesh);
			//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
			//}
		}

		TangentialRelaxation(func, mesh);
		//		AreaWeigthedTangentialRelaxation(func, mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxation.obj";
		//	obj_save(ss.str(), mesh);
		//	ss.str("");
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5_it" << i << "_AreaWeigthedTangentialRelaxations_badnormals.vtk";
		//	badNormals = VisualiseBadFaceNormals(mesh, ss.str(), true);
		//	std::cout << "Badnormals after TangentialRelaxation " << badNormals << std::endl;
		//}

		BoundarySmoothing(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step5a_it" << i << "_BoundarySmoothing.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//	std::cout << "Badnormals after BoundarySmoothing " << badNormals << std::endl;
		//}

		MaximiseMinAngle(mesh);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step6_it" << i << "_MaximiseMinAngle.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//	std::cout << "Badnormals after MaximiseMinAngle " << badNormals << std::endl;

		//}

		ProjectToSurface(func, mesh, targetLength);
		//if (debug)
		//{
		//	std::ostringstream ss;
		//	ss << debugDir + "RemeshDebug_" << saveID++ << "_Step7_it" << i << "_ProjectToSurface.obj";
		//	obj_save(ss.str(), mesh);
		//	badNormals = VisualiseBadFaceNormals(mesh, std::string(""), false);
		//}

		//std::cout << "BadNorms " << VisualiseBadFaceNormals(mesh, std::string(""), false);
		std::cout << std::endl;
	}
	std::cout << std::endl;

	locator->Delete();

	mesh.cleanup();

	//if (debug)
	//{
	//	VisualiseFaceNormals(mesh, debugDir + "RemeshDebug_Remeshed_FaceNormals.vtk");
	//	int badNormals = VisualiseBadFaceNormals(mesh, debugDir + "RemeshDebug_Remeshed_BadFaceNormals.vtk", true);
	//	std::cout << "Number of bad normals " << badNormals << std::endl;
	//}

	avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;

	GEL2VTK(mesh, pdout);
	return true;
}


bool CGELRemeshing::RemeshSphere( vtkPolyData *pdin, vtkPolyData *pdout, double radius)
{
	HMesh::Manifold mesh;

	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	if (valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
	}

	std::cout << "Computing average length" << std::endl;
	double avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;

	//	double targetLength = avgL * 0.8;
	double targetLength = avgL;
	double loTarget = 4.0/5.0   * targetLength;
	double hiTarget = 4.0 / 3.0 * targetLength;

	int iterations = 10;
	for (int i = 0; i < iterations; i++)
	{
		std::cout << "Iteration " << i << " ";

		mesh.cleanup();

		SplitLongEdges(mesh, hiTarget);

		CollapseShortEdges(mesh, loTarget, hiTarget);

		EqualizeValences(mesh);

		AreaWeigthedTangentialRelaxationOnSphere(mesh);

		MaximiseMinAngle(mesh);

		ProjectToSurfaceOnSphere(radius, mesh, targetLength);

		std::cout << std::endl;
	}
	std::cout << std::endl;

	mesh.cleanup();

	avgL = HMesh::average_edge_length(mesh);
	std::cout << "Average edge length: " << avgL << std::endl;

	GEL2VTK(mesh, pdout);
	return true;
}

void CGELRemeshing::DumpEdgeStatisticsToFile( vtkPolyData *pdin, const std::string& fname )
{
	bool BorderEdges = false;

	HMesh::Manifold mesh;

	if (!VTK2GELFaceByFace(pdin, mesh))
		return;

	std::ofstream ost(fname.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << fname << std::endl;
		return;
	}

	for(HMesh::HalfEdgeIDIterator he = mesh.halfedges_begin(); he != mesh.halfedges_end(); ++he)
	{
		if(boundary(mesh, *he) && !BorderEdges)
//		if (!HMesh::is_boundary(he) && !BorderEdges)
		{
			double l = HMesh::length(mesh, *he);
			ost << l << std::endl;
		}
	}
}

void CGELRemeshing::DumpMinAngleStatisticsToFile( vtkPolyData *pdin, const std::string& fname )
{
	std::ofstream ost(fname.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << fname << std::endl;
		return;
	}
	std::vector<double> MinAngles;
	CMeshMeasures::AllMinimumAngles(pdin, MinAngles);

	for (unsigned int i = 0; i < MinAngles.size(); i++)
	{
		ost << MinAngles[i] << std::endl;
	}

}

void CGELRemeshing::ComputeEdgeStatistics( vtkPolyData *pdin, std::vector<double> &lengths )
{
	bool BorderEdges = false;

	HMesh::Manifold mesh;

	if (!VTK2GELFaceByFace(pdin, mesh))
		return;

	for(HMesh::HalfEdgeIDIterator he = mesh.halfedges_begin(); he != mesh.halfedges_end(); ++he)
	{
		if(boundary(mesh, *he) && !BorderEdges)
//		if (!HMesh::is_boundary(he) && !BorderEdges)
		{
			double l = HMesh::length(mesh, *he);
			lengths.push_back(l);
		}
	}
}

double CGELRemeshing::ComputeAverageEdgeLength(vtkPolyData *pdin)
{
	HMesh::Manifold mesh;

	if (!VTK2GELFaceByFace(pdin, mesh))
		return 0;

	double avgL = HMesh::average_edge_length(mesh);

	return avgL;
}


void CGELRemeshing::DumpAreaStatisticsToFile( vtkPolyData *pdin, const std::string& fname )
{
	HMesh::Manifold mesh;

	if (!VTK2GELFaceByFace(pdin, mesh))
		return;

	std::ofstream ost(fname.c_str());
	if (!ost)
	{
		std::cerr << "Could not write to " << fname << std::endl;
		return;
	}

	for(HMesh::FaceIDIterator fi = mesh.faces_begin(); fi != mesh.faces_end(); ++fi)
	{
		double a = HMesh::area(mesh, *fi);
		ost << a << std::endl;
	}
}

void CGELRemeshing::MaximiseMinAngle( HMesh::Manifold &mesh )
{
	std::cout << "MaxMinAngle";
	//int flips = HMesh::maximize_min_angle(mesh, 0.9);
	//std::cout << "(" << flips << ")";
	HMesh::maximize_min_angle(mesh, 0.9);
}


// Test function
bool CGELRemeshing::Polygonise(vtkImageData * SDF, vtkImplicitFunction *func )
{
	double space[3];
	SDF->GetSpacing(space);

	double RealSS = space[0];

	vtkPolyData *pd = vtkPolyData::New();
	double cellSize = RealSS / 1;


	int resDim[3];
	SDF->GetDimensions(resDim);
	int MaxCells = std::max(resDim[0], std::max(resDim[1], resDim[2])) * 2;

	bool useTetra = false;
	std::cout << "\nBloomenthal with cell size " << cellSize << std::endl;

	double bounds[6];
	SDF->GetBounds(bounds);

	double start[3];
	start[0] = (bounds[1] + bounds[0]) / 2; 
	start[1] = (bounds[3] + bounds[2]) / 2;
	start[2] = (bounds[5] + bounds[4]) / 2;

	try
	{
		bounds[0] += 1 * cellSize;
		bounds[1] -= 1 * cellSize;
		bounds[2] += 1 * cellSize;
		bounds[3] -= 1 * cellSize;
		bounds[4] += 1 * cellSize;
		bounds[5] -= 1 * cellSize;

		Geometry::BloomenthalPolygonizer bloomenthal(func, cellSize, MaxCells, bounds, useTetra);
		bloomenthal.march(start[0], start[1], start[2]);

		using namespace HMesh;
		using namespace Geometry;
		using namespace CGLA;


		// Build HMesh from triangles
		std::vector<CGLA::Vec3d> verts;
		std::vector<int> faces;
		std::vector<int> indices;

		// Extract vertices
		for(int i=0;  i <bloomenthal.no_vertices(); ++i)
		{
			verts.push_back(*reinterpret_cast<Vec3d*>(&bloomenthal.get_vertex(i)));
		}

		// Extract triangles.
		for(int i=0;i<bloomenthal.no_triangles();++i)
		{
			faces.push_back(3);
			TRIANGLE f = bloomenthal.get_triangle(i);
			indices.push_back(f.v0);
			indices.push_back(f.v1);
			indices.push_back(f.v2);
		}

		// Build manifold
		Manifold m;
// 		build_manifold(m, verts.size(), &verts[0], 
// 			faces.size(), &faces[0], &indices[0]);
// 
//		m.build(verts.size(), reinterpret_cast<double*>(&verts[0]), faces.size(), &faces[0], &indices[0]);

		HMesh::build(m, verts.size(), reinterpret_cast<double*>(&verts[0]), faces.size(), &faces[0], &indices[0]);

		if (valid(m))
		{
			std::cout << "Custom mesh is not valid" << std::endl;
			return false;
		}


		bloomenthal.ExportToPolyData(pd);
	}
	catch (std::string ermsg)
	{
		std::cerr << "Bloomenthal fail " << ermsg << std::endl;
		return 0;
	}
	return 1;
}

bool CGELRemeshing::TestGELBoundaryOps( vtkPolyData *pdin )
{
	bool debug = false;
	std::string debugDir   = "E:/data/IMM/Remeshing/RemeshDebugDebug/";

	HMesh::Manifold mesh;

	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	if (debug)
	{
		obj_save(debugDir+"RemeshDebug.obj", mesh);
	}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	std::vector<CGLA::Vec3d> pos(mesh.allocated_vertices());
	std::vector<bool> changed(mesh.allocated_vertices(), false);

	// Run through vertices and check for boundary
	vtkPolyData *pdt = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
//	vtkCellArray *polys = vtkCellArray::New();

	int i = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi, i++)
	{
		if(boundary(mesh, *vi))
		{
			std::vector<CGLA::Vec3d> neighbours(2);

			// Test number of neighbours that are also on the boundary
			//			CGLA::Vec3d q(0);
			int n=0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				if(boundary(mesh, w.vertex()))
				{
					//					q += mesh.pos(w.vertex());
					if (n < 2)
					{
						neighbours[n] = mesh.pos(w.vertex());
					}
					n++;
				}
			}
			if (n != 2)
			{
				std::cout << "Edge neighbours " << n << std::endl;
			}

			// Compute angle between neighbouring edges
			if (n == 2)
			{
				CGLA::Vec3d v0 = mesh.pos(*vi);
				CGLA::Vec3d v1 = neighbours[0];
				CGLA::Vec3d v2 = neighbours[1];
				double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
				double Angle = vtkMath::DegreesFromRadians(a0);

				std::cout << "Angle " <<  Angle << std::endl;

				if (Angle > 170)
//				if (n != 2)
				{
					CGLA::Vec3d v = mesh.pos(*vi);

					double p[3];
					p[0] = v[0];
					p[1] = v[1];
					p[2] = v[2];

					int id = pts->InsertNextPoint(p);
					verts->InsertNextCell(1);
					verts->InsertCellPoint(id);

					pos[i] = (v1+v2) / 2;
					changed[i] = true;
				}
			}

			//if (n != 2)
			//{
			//	CGLA::Vec3d v = mesh.pos(*vi);

			//	double p[3];
			//	p[0] = v[0];
			//	p[1] = v[1];
			//	p[2] = v[2];

			//	int id = pts->InsertNextPoint(p);
			//	verts->InsertNextCell(1);
			//	verts->InsertCellPoint(id);
			//}
		}

	}

	//for(HMesh::FaceIDIterator f = mesh.faces_begin(); f != mesh.faces_end(); ++f)
	//{        
	//	std::vector<int> verts;

	//	for (HMesh::Walker w = mesh.walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
	//	{
	//		int vertex_idx = vmap[w.vertex()];			
	//		verts.push_back(vertex_idx);
	//	}
	//	if (verts.size() == 3)
	//	{
	//		polys->InsertNextCell(3);

	//		polys->InsertCellPoint(verts[0]);
	//		polys->InsertCellPoint(verts[1]);
	//		polys->InsertCellPoint(verts[2]);
	//	}
	//}

	pdt->SetPoints(pts);
	pts->Delete();
//	pdt->SetPolys(polys);
	//polys->Delete();
	pdt->SetVerts(verts);
	verts->Delete();

	if (debug)
	{
		std::string name = debugDir + "EdgeVerts.Vtk";
		vtkExtMisc::WritePDVTK(pdt, name);
	}
	
	pdt->Delete();

	i = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin();	vi != mesh.vertices_end(); ++vi, i++)
	{
		if (changed[i])
		{
			mesh.pos(*vi) = pos[i];
		}
	}

	if (debug)
	{
		obj_save(debugDir+"RemeshDebugBoundaryRelaxed.obj", mesh);
	}


	return true;
}


bool CGELRemeshing::TestGELBoundaryCollapseEdge( vtkPolyData *pdin )
{
	bool debug = false;
	std::string debugDir   = "E:/data/IMM/Remeshing/RemeshDebugDebug/";

	HMesh::Manifold mesh;

	std::cout << "VTK2GELFaceByFace" << std::endl;
	if (!VTK2GELFaceByFace(pdin, mesh))
		return false;

	if (debug)
	{
		obj_save(debugDir+"RemeshDebug.obj", mesh);
	}

	if (!valid(mesh))
	{
		std::cout << "GEL claims that mesh is not valid" << std::endl;
		std::cout << "Still trying to remesh" << std::endl;
		//		return false;
	}

	vtkPolyData *pdt = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();

	int i = 0;
	for(HMesh::VertexIDIterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi, i++)
	{
		if(boundary(mesh, *vi))
		{
			std::vector<CGLA::Vec3d> neighbours(2);

			int n=0;
			HMesh::Walker w = mesh.walker(*vi);
			for(; !w.full_circle(); w = w.circulate_vertex_ccw())
			{
				if(boundary(mesh, w.vertex()))
				{
					if (n < 2)
					{
						neighbours[n] = mesh.pos(w.vertex());
					}
					n++;
				}
			}

			double low = 1;

			// Compute angle between neighbouring edges
			if (n == 2)
			{
				CGLA::Vec3d v0 = mesh.pos(*vi);
				CGLA::Vec3d v1 = neighbours[0];
				CGLA::Vec3d v2 = neighbours[1];
				double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
				double Angle = vtkMath::DegreesFromRadians(a0);

				if (Angle > 170)
				{
					HMesh::HalfEdgeID h = w.halfedge();
					double dist = length(mesh, h);
					{
						// collapse edge if allowed
						if(dist < low && precond_collapse_edge(mesh, h))
						{
							std::cout << "Collapsing boundary edge of length " << dist << std::endl;
							CGLA::Vec3d v = mesh.pos(*vi);
							double p[3];
							p[0] = v[0];
							p[1] = v[1];
							p[2] = v[2];

							int id = pts->InsertNextPoint(p);
							verts->InsertNextCell(1);
							verts->InsertCellPoint(id);

							mesh.collapse_edge(h);
						}
					}
				}
			}
		}

	}

	pdt->SetPoints(pts);
	pts->Delete();
	pdt->SetVerts(verts);
	verts->Delete();

	if (debug)
	{
		std::string name = debugDir + "ShortEdgeVerts.Vtk";
		vtkExtMisc::WritePDVTK(pdt, name);
	}

	pdt->Delete();

	if (debug)
	{
		obj_save(debugDir+"RemeshDebugBoundaryCollapsed.obj", mesh);
	}

	return true;
}




