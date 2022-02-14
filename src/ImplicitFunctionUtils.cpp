#include "ImplicitFunctionUtils.h"

#include <vtkImplicitFunction.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <math.h>

CImplicitFunctionUtils::CImplicitFunctionUtils()
{
}

CImplicitFunctionUtils::~CImplicitFunctionUtils()
{
}

bool CImplicitFunctionUtils::FindClosestPointByBisection(vtkImplicitFunction* func, double* p, double* cp, double epsilon)
{
	double* n = func->FunctionGradient(p);
	double dist = func->FunctionValue(p);

	// Check gradient and distance
	if (n[0] == 0 && n[1] == 0 && n[2] == 0)
	{
		return false;
	}
	if (std::abs(dist) > 1e100)
	{
		return false;
	}

	double a = 0;
	double b = dist * 3;
	double minPos = a;
	double fx1 = epsilon * 10;
	double fx2 = epsilon * 10;
	int its = 0;
	do
	{
		double x1 = (b - a) / 3 + a;
		double x2 = 2 * (b - a) / 3 + a;

		fx1 = func->FunctionValue(p[0] + x1 * n[0], p[1] + x1 * n[1], p[2] + x1 * n[2]);
		fx2 = func->FunctionValue(p[0] + x2 * n[0], p[1] + x2 * n[1], p[2] + x2 * n[2]);

		if (fx1 > fx2)
		{
			a = x1;
			minPos = x2;
		}
		else
		{
			b = x2;
			minPos = x1;
		}
	} while (std::min(fx1, fx2) > epsilon && its++ < 20);

	cp[0] = p[0] + minPos * n[0];
	cp[1] = p[1] + minPos * n[1];
	cp[2] = p[2] + minPos * n[2];

	if (std::abs(cp[0]) > 1000000 || std::abs(cp[1]) > 1000000 || std::abs(cp[2]) > 1000000)
	{
		std::cout << "Point p(" << p[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;
		std::cout << "Closest Point cp(" << cp[0] << ", " << cp[1] << ", " << cp[2] << ")" << std::endl;
	}

	return true;
}

bool CImplicitFunctionUtils::FindClosestPoint( vtkImplicitFunction *func, double *p, double *cp, double accuracy, double SanityJump, double MaxInitialDist /*= 1e8*/, int iterations)
{
	double dist;
	double *n;
	cp[0] = p[0]; cp[1] = p[1]; cp[2] = p[2];

	double thres = accuracy;
	double fac = 0.5;
	//double fac = 0.5;

	dist = func->FunctionValue(cp);
	if (abs(dist) > MaxInitialDist)
	{
		//std::cerr << "Not ok in ProjectToSurface:abs(dist)" << abs(dist) << " > MaxInitialDist " << MaxInitialDist << std::endl;
		return false;
	}
	
	if (abs(dist) <= thres)
	{
	//	std::cerr << "Ok in ProjectToSurface: abs(dist) <= thres" << std::endl;
		return true;
	}


	n = func->FunctionGradient(cp);
	double nl = sqrt(n[0] * n[0]  + n[1] * n[1] + n[2] * n[2]);
	if (nl)
	{
		n[0] /= nl;
		n[1] /= nl;
		n[2] /= nl;
	}

//	std::cout << "Initial distance " << dist << std::endl;

	double oldDist = abs(dist);

	double minDist = oldDist;
	double minP[3];
	minP[0] = cp[0];
	minP[1] = cp[1];
	minP[2] = cp[2];
	int NoMinIt = 0;

	int nsteps = 0;
	do 
	{
		cp[0] = cp[0] + n[0] * fac * dist;
		cp[1] = cp[1] + n[1] * fac * dist;
		cp[2] = cp[2] + n[2] * fac * dist;

		dist = func->FunctionValue(cp);

		// This typically happens at the edges of the volume
		if (abs(dist) > oldDist * 2)
		{
			//std::cerr << "Not ok in ProjectToSurface: abs(dist): " << abs(dist) << " > oldDist * 2: " << oldDist * 2 <<
			//	 " for point : (" << cp[0] << ", " << cp[1] << ", " << cp[2] << ") dist: " << dist << 
			//	 " n : (" << n[0] << ", " << n[1] << ", " << n[2] << std::endl;

			return false;
		}
		if (abs(dist) > 1e8)
		{
			//std::cerr << "Not ok in ProjectToSurface: abs(dist) > 1e8" << std::endl;
//			std::cerr << "Something wrong in distance. Dist " << dist << std::endl;
			return false;
		}
		if (abs(dist) < minDist)
		{
			minDist = abs(dist);
			minP[0] = cp[0];
			minP[1] = cp[1];
			minP[2] = cp[2];
			NoMinIt = 0;
		}
		else if (NoMinIt++ > iterations / 2)
		{
			//if (minDist > 2)
			//	std::cout << "Oscillated for more than : " << iterations / 2 << " nsteps: " << nsteps << " dist " << dist << " olddist " << oldDist << " mindist " << minDist << " nl " << nl << " cp (" << cp[0] << ", " << cp[1] << ", " << cp[2] << ")" << std::endl;

			// No new minimum found in 5 iterations
			cp[0] = minP[0];
			cp[1] = minP[1];
			cp[2] = minP[2];
			return true;
		}

		n = func->FunctionGradient(cp);
		double nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (nl)
		{
			n[0] /= nl;
			n[1] /= nl;
			n[2] /= nl;
		}

		nsteps++;
		//if (nsteps > 500)
		//{
		//	std::cout << "ProjectToSurface: " << nsteps << " dist " << dist << " nl " << nl << " cp (" << cp[0] << ", " << cp[1] << ", " << cp[2] << ")" << std::endl;
		//}

		if (nsteps > iterations)
		{
			//std::cerr << "Not ok in ProjectToSurface: nsteps > 1000" << std::endl;
			//if (minDist > 2)
			//	std::cout << "ProjectToSurface: " << nsteps << " dist " << dist << " olddist " << oldDist << " mindist " << minDist << " nl " << nl << " cp (" << cp[0] << ", " << cp[1] << ", " << cp[2] << ")" << std::endl;

			cp[0] = minP[0];
			cp[1] = minP[1];
			cp[2] = minP[2];
			return true;
		}

		oldDist = abs(dist);

//		std::cout << "Distance " << dist << std::endl;
	} while (abs(dist) > thres);
	// std::cout << "dist: " << dist << std::endl;

	// Check jump distance. Probably mostly happens at the edge of the volume
	double JumpDist = (p[0]-cp[0])*(p[0]-cp[0]) + (p[1]-cp[1])*(p[1]-cp[1]) + (p[2]-cp[2])*(p[2]-cp[2]);
	if (JumpDist  > SanityJump * SanityJump)
	{
		//std::cerr << "Not ok in ProjectToSurface: Jump: " << JumpDist << " longer than sanity jump. dist " << dist << " mindist " << minDist << std::endl;
		return false;
	}

	return true;
}

bool CImplicitFunctionUtils::FindClosestPointIllustrated( vtkImplicitFunction *func, double *p, double *cp, vtkPolyData *path, double accuracy)
{
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	double dist;
	double *n;
	cp[0] = p[0]; cp[1] = p[1]; cp[2] = p[2];

	double thres = 0.1;
	double fac = 0.5;

	dist = func->FunctionValue(cp);
	n = func->FunctionGradient(cp);

	std::cout << "Initial distance " << dist << std::endl;

	verts->InsertNextCell(1);
	vtkIdType id = pts->InsertNextPoint(cp);
	verts->InsertCellPoint(id);
	scalars->InsertNextTuple(&dist);


	do 
	{
		cp[0] = cp[0] + n[0] * fac * dist;
		cp[1] = cp[1] + n[1] * fac * dist;
		cp[2] = cp[2] + n[2] * fac * dist;

		dist = func->FunctionValue(cp);
		n = func->FunctionGradient(cp);

		verts->InsertNextCell(1);
		vtkIdType id = pts->InsertNextPoint(cp);
		verts->InsertCellPoint(id);
		scalars->InsertNextTuple(&dist);

		std::cout << "Distance " << dist << std::endl;
	} while (abs(dist) > thres);

	path->SetPoints(pts);
	pts->Delete();
	path->SetVerts(verts);
	verts->Delete();
	path->GetPointData()->SetScalars(scalars);
	scalars->Delete();

	return true;
}
