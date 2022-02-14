#include "MeshMeasures.h"

#include <vtkpolydata.h>
#include <vtkPointLocator.h>
#include "vtkMath.h"
#include "vtkCell.h"
#include "vtkIdList.h"
#include "vtkTriangle.h"

#include <iostream>
#include "GeneralUtils.h"

void CMeshMeasures::MeshRegularityBasedOnEdgeLengths(vtkPolyData *mesh, double &DistMean, double &DistSdev)
{
	mesh->BuildLinks();


	int NP = mesh->GetNumberOfPoints();

	double S = 0.0;
	double SS = 0.0;
	int N = 0;
	for (int i = 0; i < NP; i++)
	{
		for (int j = i+1; j < NP; j++)
		{
			if (mesh->IsEdge(i,j))
			{
				double pNi[3]; double pNj[3];
				mesh->GetPoint(i, pNi); mesh->GetPoint(j, pNj);
				double d = sqrt(vtkMath::Distance2BetweenPoints(pNi, pNj));
				S  += d;
				SS += d*d;
				N++;
			}
		}
	}

	DistMean = 0;
	DistSdev = 0;
	if (N)
	{
		DistMean = S/N;
		DistSdev = sqrt(SS/N - DistMean*DistMean);
	}
}

void CMeshMeasures::MeshRegularityBasedOnMinimumTriangleAngle(vtkPolyData *mesh, double &MinAngleMean)
{
	int NP = mesh->GetNumberOfCells();

	double S = 0.0;
	int N = 0;
	for (int i = 0; i < NP; i++)
	{
		vtkTriangle *tri = vtkTriangle::SafeDownCast(mesh->GetCell(i));
		if (!tri)
		{
			std::cerr << "Cell " << i << " is not a triangle!" << std::endl;
		}
		else
		{
			
			int id0 = tri->GetPointId(0);
			int id1 = tri->GetPointId(1);
			int id2 = tri->GetPointId(2);

			// cos A = (b*b + c*c - a*a) / (2 * b *c)

			double pID0[3]; double pID1[3]; double pID2[3];
			mesh->GetPoint(id0, pID0); mesh->GetPoint(id1, pID1); mesh->GetPoint(id2, pID2);

			double a2 = vtkMath::Distance2BetweenPoints(pID0, pID1);
			double b2 = vtkMath::Distance2BetweenPoints(pID1, pID2);
			double c2 = vtkMath::Distance2BetweenPoints(pID0, pID2);

			if (a2 == 0 || b2 == 0 || c2 == 0)
			{
				std::cerr << "Triangle " << i << " degenerate. Sidelength 0" << std::endl;
			}
			else
			{
				double a = sqrt(a2);
				double b = sqrt(b2);
				double c = sqrt(c2);

				double cosA = (b2 + c2 - a2) / (2 * b * c);
				double cosB = (c2 + a2 - b2) / (2 * c * a);
				double cosC = (a2 + b2 - c2) / (2 * a * b);

				double A = vtkMath::DegreesFromRadians(fabs(acos(cosA)));
				double B = vtkMath::DegreesFromRadians(fabs(acos(cosB)));
				double C = vtkMath::DegreesFromRadians(fabs(acos(cosC)));
		
				double totalA = A + B + C;

				double minA = std::min(A, std::min(B,C));

				S += minA;
				N++;
			}
		}
	}
	if (N)
	{
		S /= N;
	}
	MinAngleMean = S;
}


bool CMeshMeasures::IsCellContainedInOtherCell(vtkCell *cell1, vtkCell *cell2)
{
	vtkIdList *ids1 = cell1->GetPointIds();
	vtkIdList *ids2 = cell2->GetPointIds();

	// Check if all ids in ids1 is also in ids2
	for (int i = 0; i < ids1->GetNumberOfIds(); i++)
	{
		vtkIdType id1 = ids1->GetId(i);

		bool found = false;
		for (int j = 0; j < ids2->GetNumberOfIds() && !found; j++)
		{
			if (id1 == ids2->GetId(j))
			{
				found = true;
			}
		}
		if (found == false)
		{
			return false;
		}
	}

	return true;
}

bool CMeshMeasures::IsCellContainedInOtherCell(vtkIdType npts1, const vtkIdType *pts1, vtkIdType npts2, const vtkIdType *pts2)
{
	// Check if all ids in ids1 is also in ids2
	for (int i = 0; i < npts1; i++)
	{
		vtkIdType id1 = pts1[i];

		bool found = false;
		for (int j = 0; j < npts2 && !found; j++)
		{
			if (id1 == pts2[j])
			{
				found = true;
			}
		}
		if (found == false)
		{
			return false;
		}
	}

	return true;
}


int CMeshMeasures::TotalNumberOfCellsContainedInOtherCells(vtkPolyData *pd)
{
	pd->BuildCells();
	const int NCells = pd->GetNumberOfCells();

	int count = 0;
	for (int i = 0; i < NCells; i++)
	{
		if (!(i % 500)) std::cout << i << std::endl;

		for (int j = 0; j < NCells; j++)
		{
			if (i != j)
			{
				vtkIdType npts1 = 0;
				vtkIdType npts2 = 0;
				const vtkIdType *pts1 = NULL;
				const vtkIdType *pts2 = NULL;
		
				pd->GetCellPoints(i, npts1, pts1);
				pd->GetCellPoints(j, npts2, pts2);

				if (CMeshMeasures::IsCellContainedInOtherCell(npts1, pts1, npts2, pts2))
				{
					count++;
				}
			}
		}
	}
	return count;
}

void CMeshMeasures::AllMinimumAngles( vtkPolyData *mesh, std::vector<double> &MinAngles )
{
	MinAngles.clear();

	int NP = mesh->GetNumberOfCells();

	for (int i = 0; i < NP; i++)
	{
		vtkTriangle *tri = vtkTriangle::SafeDownCast(mesh->GetCell(i));
		if (!tri)
		{
			std::cerr << "Cell " << i << " is not a triangle!" << std::endl;
		}
		else
		{

			int id0 = tri->GetPointId(0);
			int id1 = tri->GetPointId(1);
			int id2 = tri->GetPointId(2);

			// cos A = (b*b + c*c - a*a) / (2 * b *c)

			double pID0[3]; double pID1[3]; double pID2[3];
			mesh->GetPoint(id0, pID0); mesh->GetPoint(id1, pID1); mesh->GetPoint(id2, pID2);

			double a2 = vtkMath::Distance2BetweenPoints(pID0, pID1);
			double b2 = vtkMath::Distance2BetweenPoints(pID1, pID2);
			double c2 = vtkMath::Distance2BetweenPoints(pID0, pID2);

			if (a2 == 0 || b2 == 0 || c2 == 0)
			{
				std::cerr << "Triangle " << i << " degenerate. Sidelength 0" << std::endl;
			}
			else
			{
				double a = sqrt(a2);
				double b = sqrt(b2);
				double c = sqrt(c2);

				double cosA = (b2 + c2 - a2) / (2 * b * c);
				double cosB = (c2 + a2 - b2) / (2 * c * a);
				double cosC = (a2 + b2 - c2) / (2 * a * b);

				double A = vtkMath::DegreesFromRadians(fabs(acos(cosA)));
				double B = vtkMath::DegreesFromRadians(fabs(acos(cosB)));
				double C = vtkMath::DegreesFromRadians(fabs(acos(cosC)));

	//			double totalA = A + B + C;

				double minA = std::min(A, std::min(B,C));

				MinAngles.push_back(minA);
			}
		}
	}
}


bool CMeshMeasures::NearestNeighbourStatistics(vtkDataSet *pd, double &minDist, double &maxDist, double &meanDist, double &sdvDist, double &medianDist, double &frac05, double &frac95)
{
	const int NPoint = pd->GetNumberOfPoints();
	if (NPoint < 2)
	{
		std::cerr << "Less than 2 points" << std::endl;
		return false;
	}

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(pd);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	std::vector<double> Tdists;

	for (int i = 0; i < NPoint; i++)
	{
		double p[3];
		pd->GetPoint(i, p);

		vtkIdList *neighPts = vtkIdList::New();
		locator->FindClosestNPoints(2, p, neighPts);

		if (neighPts->GetNumberOfIds() == 2)
		{
			vtkIdType cid = neighPts->GetId(1);
			double cp[3];

			pd->GetPoint(cid, cp);


			double dist = sqrt(vtkMath::Distance2BetweenPoints(p, cp));
			if (dist > 0)
			{
				Tdists.push_back(dist);
			}
		}

		neighPts->Delete();
	}
	locator->Delete();

	
	CGeneralUtils::MeanAndSdev(Tdists, meanDist, sdvDist);
	CGeneralUtils::MinMax(Tdists, minDist, maxDist);
	CGeneralUtils::Median(Tdists, 0.05, frac05);
	CGeneralUtils::Median(Tdists, 0.95, frac95);
	CGeneralUtils::Median(Tdists, 0.50, medianDist);
	return true;
}

bool CMeshMeasures::TriangleSideLengths(vtkPolyData *mesh, int cellId, double &l1, double &l2, double &l3)
{
	vtkTriangle *tri = vtkTriangle::SafeDownCast(mesh->GetCell(cellId));
	if (!tri)
	{
		std::cerr << "Cell " << cellId << " is not a triangle!" << std::endl;
		return false;
	}

	int id0 = tri->GetPointId(0);
	int id1 = tri->GetPointId(1);
	int id2 = tri->GetPointId(2);

	double pID0[3]; double pID1[3]; double pID2[3];
	mesh->GetPoint(id0, pID0); mesh->GetPoint(id1, pID1); mesh->GetPoint(id2, pID2);

	double a2 = vtkMath::Distance2BetweenPoints(pID0, pID1);
	double b2 = vtkMath::Distance2BetweenPoints(pID1, pID2);
	double c2 = vtkMath::Distance2BetweenPoints(pID0, pID2);

	if (a2 == 0 || b2 == 0 || c2 == 0)
	{
//		std::cerr << "Triangle " << i << " degenerate. Sidelength 0" << std::endl;
		return false;
	}
	
	l1 = sqrt(a2);
	l2 = sqrt(b2);
	l3 = sqrt(c2);
	return true;
}

bool CMeshMeasures::TriangleAngles( vtkPolyData *mesh, int cellId, double &A, double &B, double &C, double &minA )
{
	vtkTriangle *tri = vtkTriangle::SafeDownCast(mesh->GetCell(cellId));
	if (!tri)
	{
		std::cerr << "Cell " << cellId << " is not a triangle!" << std::endl;
		return false;
	}

	int id0 = tri->GetPointId(0);
	int id1 = tri->GetPointId(1);
	int id2 = tri->GetPointId(2);

	double pID0[3]; double pID1[3]; double pID2[3];
	mesh->GetPoint(id0, pID0); mesh->GetPoint(id1, pID1); mesh->GetPoint(id2, pID2);

	double a2 = vtkMath::Distance2BetweenPoints(pID0, pID1);
	double b2 = vtkMath::Distance2BetweenPoints(pID1, pID2);
	double c2 = vtkMath::Distance2BetweenPoints(pID0, pID2);

	if (a2 == 0 || b2 == 0 || c2 == 0)
	{
		//		std::cerr << "Triangle " << i << " degenerate. Sidelength 0" << std::endl;
		return false;
	}

	double a = sqrt(a2);
	double b = sqrt(b2);
	double c = sqrt(c2);

	double cosA = (b2 + c2 - a2) / (2 * b * c);
	double cosB = (c2 + a2 - b2) / (2 * c * a);
	double cosC = (a2 + b2 - c2) / (2 * a * b);

	A = vtkMath::DegreesFromRadians(fabs(acos(cosA)));
	B = vtkMath::DegreesFromRadians(fabs(acos(cosB)));
	C = vtkMath::DegreesFromRadians(fabs(acos(cosC)));

	minA = std::min(A, std::min(B,C));
	return true;
}
