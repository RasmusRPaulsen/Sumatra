#include "vtkPolyDataPCANormals.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
//#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
//#include "vtkPolygon.h"
//#include "vtkTriangleStrip.h"
//#include "vtkPriorityQueue.h"

#include <vtkTimerLog.h>

#include <vtkDoubleArray.h>
#include <vtkPointLocator.h>

vtkStandardNewMacro(vtkPolyDataPCANormals);


vtkPolyDataPCANormals::vtkPolyDataPCANormals()
{
	MaxPlaneDistance = 1.0;
	SearchRadius = 1.0;
	Max3EigenValue = 20.0;
//	ScalarMode = VTK_SCALAR_MODE_3RDEIGENVALUE;
	ScalarMode = VTK_SCALAR_MODE_PLANEDIST;
	PointsPerNormals = 5;
	SearchMode = VTK_SEARCH_MODE_SEARCHRADIUS;
}

void vtkPolyDataPCANormals::CreateConnectivityByPointsPerNormal(vtkPolyData *input)
{
	std::cout << "Creating connectivity by searching for " << PointsPerNormals << " points per normal" << std::endl;

	m_Connectivity.resize(input->GetNumberOfPoints());

	vtkPointLocator *Locator = vtkPointLocator::New();
	Locator->SetDataSet(input);
	Locator->SetNumberOfPointsPerBucket(1);
	Locator->BuildLocator();

	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double p[3];
		input->GetPoint(i, p);

		vtkIdList *idlist = vtkIdList::New();
		Locator->FindClosestNPoints(PointsPerNormals, p, idlist);

		for (int j = 0; j < idlist->GetNumberOfIds(); j++)
		{
			int idx = idlist->GetId(j);
			if (i != idx)
			{
				m_Connectivity[i].push_back(idx);
			}
		}

		idlist->Delete();
	}

	Locator->Delete();
}

void vtkPolyDataPCANormals::CreateConnectivityUsingSearchRadius(vtkPolyData *input)
{
	std::cout << "Creating connectivity by searching in a radius of " << SearchRadius << "  per normal" << std::endl;

	m_Connectivity.resize(input->GetNumberOfPoints());

	vtkPointLocator *Locator = vtkPointLocator::New();
	Locator->SetDataSet(input);
	Locator->SetNumberOfPointsPerBucket(1);
	Locator->BuildLocator();


	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double p[3];
		input->GetPoint(i, p);

		vtkIdList *idlist = vtkIdList::New();
		Locator->FindPointsWithinRadius(SearchRadius, p, idlist);

		for (int j = 0; j < idlist->GetNumberOfIds(); j++)
		{
			int idx = idlist->GetId(j);
			if (i != idx)
			{
				m_Connectivity[i].push_back(idx);
			}
		}

		idlist->Delete();
	}

	Locator->Delete();
}

void ComputeEigenVectorsOfPointCloud(vtkPolyData *pd, std::vector<int>& ids, double* mean, 
												 double *v0, double *v1, double *v2, double* evals)
{
	unsigned int i,j;

//	vtkMatrix4x4 *Matrix = vtkMatrix4x4::New();

	if (ids.size() == 0)
	{
		std:: cerr << "Can't execute with NULL or empty input";
		return;
	}

	size_t numPts = ids.size();
	mean[0] = 0; mean[1] = 0; mean[2] = 0;
	for (i = 0; i < ids.size(); i++)
	{
		int idx = ids[i];

		double x[3];
		pd->GetPoint(idx, x);
		mean[0] += x[0];
		mean[1] += x[1];
		mean[2] += x[2];
	}
	mean[0] /= numPts;
	mean[1] /= numPts;
	mean[2] /= numPts;


	// a is a vector of pointers pointing into the three vectors a0, a1 and a2
	// so basically a is a matrix with rows a0, a1 and a2
	double *a[3], a0[3], a1[3], a2[3];

	a[0] = a0; a[1] = a1; a[2] = a2; 
	for (i=0; i < 3; i++)
	{
		a0[i] = a1[i] = a2[i] = 0.0;
	}

	// Compute covariance matrix
	//
	for (j = 0; j < ids.size(); j++)
	{
		vtkIdType pointId = ids[j];
	
		double x[3];
		pd->GetPoint(pointId, x);

		double xp[3];
		xp[0] = x[0] - mean[0]; xp[1] = x[1] - mean[1]; xp[2] = x[2] - mean[2];

		for (i=0; i < 3; i++)
		{
			a0[i] += xp[0] * xp[i];
			a1[i] += xp[1] * xp[i];
			a2[i] += xp[2] * xp[i];
		}
	}//for all points

	for (i=0; i < 3; i++)
	{
		a0[i] /= numPts;
		a1[i] /= numPts;
		a2[i] /= numPts;
	}

	//
	// Extract axes (i.e., eigenvectors) from covariance matrix. 
	//
	double *v[3]; //, v0[3], v1[3], v2[3];
	//	double evals[3] = {0.0, 0.0, 0.0};
	evals[0] = 0; evals[1] = 0; evals[2] = 0;
	// v is a vector of pointers pointing into the three vectors v0, v1 and v2
	// so basically v is a matrix with rows v0, v1 and v2
	// the eigenvectors are returned in v. the eigenvalues in evals;

	v[0] = v0; v[1] = v1; v[2] = v2; 
	vtkMath::Jacobi(a,evals,v);
}

bool ComputeEigenVectorsOfPointCloud(vtkPolyData *pd, vtkIdList *ids, double* mean, 
									 double *v0, double *v1, double *v2, double* evals, vtkIdType id)
{
	int i,j;

	if (ids->GetNumberOfIds() == 0)
	{
		std:: cerr << "Can't execute with NULL or empty input";
		return false;
	}

	int numPts = ids->GetNumberOfIds();
	mean[0] = 0; mean[1] = 0; mean[2] = 0;
	for (i = 0; i < ids->GetNumberOfIds(); i++)
	{
		int idx = ids->GetId(i);
		if (idx == id)
		{
			numPts--;
		}
		else
		{
			double x[3];
			pd->GetPoint(idx, x);
			mean[0] += x[0];
			mean[1] += x[1];
			mean[2] += x[2];
		}
	}
	mean[0] /= numPts;
	mean[1] /= numPts;
	mean[2] /= numPts;


	// a is a vector of pointers pointing into the three vectors a0, a1 and a2
	// so basically a is a matrix with rows a0, a1 and a2
	double *a[3], a0[3], a1[3], a2[3];

	a[0] = a0; a[1] = a1; a[2] = a2; 
	for (i=0; i < 3; i++)
	{
		a0[i] = a1[i] = a2[i] = 0.0;
	}

	// Compute covariance matrix
	//
	for (j = 0; j < ids->GetNumberOfIds(); j++)
	{
		vtkIdType pointId = ids->GetId(j);

		if (pointId != id)
		{
			double x[3];
			pd->GetPoint(pointId, x);

			double xp[3];
			xp[0] = x[0] - mean[0]; xp[1] = x[1] - mean[1]; xp[2] = x[2] - mean[2];

			for (i=0; i < 3; i++)
			{
				a0[i] += xp[0] * xp[i];
				a1[i] += xp[1] * xp[i];
				a2[i] += xp[2] * xp[i];
			}
		}

	}//for all points

	for (i=0; i < 3; i++)
	{
		a0[i] /= numPts;
		a1[i] /= numPts;
		a2[i] /= numPts;
	}

	//
	// Extract axes (i.e., eigenvectors) from covariance matrix. 
	//
	double *v[3]; //, v0[3], v1[3], v2[3];
	//	double evals[3] = {0.0, 0.0, 0.0};
	evals[0] = 0; evals[1] = 0; evals[2] = 0;
	// v is a vector of pointers pointing into the three vectors v0, v1 and v2
	// so basically v is a matrix with rows v0, v1 and v2
	// the eigenvectors are returned in v. the eigenvalues in evals;

	v[0] = v0; v[1] = v1; v[2] = v2; 
	vtkMath::Jacobi(a,evals,v);

	return true;
}

void vtkPolyDataPCANormals::ComputeLocalPCAByPointsPerNormal( vtkPolyData *input, vtkPolyData *output )
{
	std::cout << "Compute local PCA by searching in a radius of " << SearchRadius << "  per normal" << std::endl;

	vtkPointLocator *Locator = vtkPointLocator::New();
	Locator->SetDataSet(input);
	Locator->SetNumberOfPointsPerBucket(1);
	Locator->BuildLocator();

	int NNoNorms = 0;
	int GoodNorms = 0;
	int badEvals = 0;
	int noisePoints = 0;
	int distPoints = 0;
	int MinPointUsed = 1000000;
	int MaxPointUsed = 0;
	int TotalPointUsed = 0;

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *Normals = vtkDoubleArray::New();
	Normals->SetNumberOfComponents(3);
	Normals->SetName("Normals");
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double mean[3];
		double v0[3];
		double v1[3];
		double v2[3];
		double evals[3];

		double p[3];
		input->GetPoint(i, p);

		vtkIdList *idlist = vtkIdList::New();
		Locator->FindClosestNPoints(PointsPerNormals, p, idlist);

		if (idlist->GetNumberOfIds() > 3)
		{
			bool result = ComputeEigenVectorsOfPointCloud(input, idlist, mean, v0, v1, v2, evals, i);

			if ((evals[0] + evals[1] + evals[2]) > 0)
			{
				// variance explained by 3rd eigenvector
				double varexp = varexp = evals[2] / (evals[0] + evals[1] + evals[2]) * 100.0;

				if (varexp < Max3EigenValue)
				{
					double n[3];
					n[0] = v0[2];
					n[1] = v1[2];
					n[2] = v2[2];
					vtkMath::Normalize(n);

					double p[3];
					input->GetPoint(i, p);

					double pt[3];
					pt[0] = p[0] - mean[0];
					pt[1] = p[1] - mean[1];
					pt[2] = p[2] - mean[2];

					double planeDist = abs(vtkMath::Dot(pt, n));

					if (planeDist < MaxPlaneDistance)
					{
						vtkIdType id = pts->InsertNextPoint(p);
						verts->InsertNextCell(1);
						verts->InsertCellPoint(id);

						Normals->InsertNextTuple(n);

						if (ScalarMode == VTK_SCALAR_MODE_3RDEIGENVALUE)
						{
							scalars->InsertNextTuple(&varexp);
						}
						else // distance to plane
						{
							scalars->InsertNextTuple(&planeDist);
						}
						GoodNorms++;
					}
					else
					{
						distPoints++;
					}
				}
				else
				{
					noisePoints++;
				}
			}
			else
			{
				badEvals++;
			}
		}
		else
		{
			NNoNorms++;
		}
		idlist->Delete();
	}

	output->SetPoints(pts);
	output->SetVerts(verts);
	output->GetPointData()->SetScalars(scalars);
	output->GetPointData()->SetNormals(Normals);

	Normals->Delete();
	pts->Delete();
	verts->Delete();
	scalars->Delete();
	Locator->Delete();

	std::cout << "Created " << GoodNorms << " points with normals" << std::endl;
	std::cout << "Points with not enough neighbours to use for normal calculation " << NNoNorms << std::endl;
	std::cout << "Points with bad evals " <<  badEvals << std::endl;
	std::cout << "Points where the 3rd Eigenvector explained more than " << Max3EigenValue << " % of the variation " << noisePoints << std::endl;
	std::cout << "Points where the distance to the fitted plane was more than " << MaxPlaneDistance << " mm " << distPoints << std::endl;
}


void vtkPolyDataPCANormals::ComputeLocalPCAUsingSearchRadius( vtkPolyData *input, vtkPolyData *output )
{
	std::cout << "Compute local PCA by searching in a radius of " << SearchRadius << "  per normal" << std::endl;

	vtkPointLocator *Locator = vtkPointLocator::New();
	Locator->SetDataSet(input);
	Locator->SetNumberOfPointsPerBucket(1);
	Locator->BuildLocator();

	int NNoNorms = 0;
	int GoodNorms = 0;
	int badEvals = 0;
	int noisePoints = 0;
	int distPoints = 0;
	int MinPointUsed = 1000000;
	int MaxPointUsed = 0;
	int TotalPointUsed = 0;

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *Normals = vtkDoubleArray::New();
	Normals->SetNumberOfComponents(3);
	Normals->SetName("Normals");
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double mean[3];
		double v0[3];
		double v1[3];
		double v2[3];
		double evals[3];
		double p[3];
		
		input->GetPoint(i, p);

		vtkIdList *idlist = vtkIdList::New();
		Locator->FindPointsWithinRadius(SearchRadius, p, idlist);

		if (idlist->GetNumberOfIds() > 3)
		{
			bool result = ComputeEigenVectorsOfPointCloud(input, idlist, mean, v0, v1, v2, evals, i);

			if ((evals[0] + evals[1] + evals[2]) > 0)
			{
				// variance explained by 3rd eigenvector
				double varexp = varexp = evals[2] / (evals[0] + evals[1] + evals[2]) * 100.0;

				if (varexp < Max3EigenValue)
				{
					double n[3];
					n[0] = v0[2];
					n[1] = v1[2];
					n[2] = v2[2];
					vtkMath::Normalize(n);

					double p[3];
					input->GetPoint(i, p);

					double pt[3];
					pt[0] = p[0] - mean[0];
					pt[1] = p[1] - mean[1];
					pt[2] = p[2] - mean[2];

					double planeDist = abs(vtkMath::Dot(pt, n));

					if (planeDist < MaxPlaneDistance)
					{
						vtkIdType id = pts->InsertNextPoint(p);
						verts->InsertNextCell(1);
						verts->InsertCellPoint(id);

						Normals->InsertNextTuple(n);

						if (ScalarMode == VTK_SCALAR_MODE_3RDEIGENVALUE)
						{
							scalars->InsertNextTuple(&varexp);
						}
						else // distance to plane
						{
							scalars->InsertNextTuple(&planeDist);
						}
						GoodNorms++;

						int realNumPoints = idlist->GetNumberOfIds() - 1; // not using current point
						if (realNumPoints< MinPointUsed)
							MinPointUsed = realNumPoints;
						if (realNumPoints > MaxPointUsed)
							MaxPointUsed = realNumPoints;
						TotalPointUsed += realNumPoints;
					}
					else
					{
						distPoints++;
					}
				}
				else
				{
					noisePoints++;
				}
			}
			else
			{
				badEvals++;
			}
		}
		else
		{
			NNoNorms++;
		}
		idlist->Delete();
	}

	output->SetPoints(pts);
	output->SetVerts(verts);
	output->GetPointData()->SetScalars(scalars);
	output->GetPointData()->SetNormals(Normals);

	Normals->Delete();
	pts->Delete();
	verts->Delete();
	scalars->Delete();
	Locator->Delete();

	std::cout << "Created " << GoodNorms << " points with normals" << std::endl;
	std::cout << "Points with not enough neighbours to use for normal calculation " << NNoNorms << std::endl;
	std::cout << "Points with bad evals " <<  badEvals << std::endl;
	std::cout << "Points where the 3rd Eigenvector explained more than " << Max3EigenValue << " % of the variation " << noisePoints << std::endl;
	std::cout << "Points where the distance to the fitted plane was more than " << MaxPlaneDistance << " mm " << distPoints << std::endl;
	if (GoodNorms)
		std::cout << "For the normals computed: min number of point used: " << MinPointUsed << " max number " << MaxPointUsed << " average number " << TotalPointUsed / GoodNorms << std::endl;
}

void vtkPolyDataPCANormals::ComputeLocalPCA(vtkPolyData *input, vtkPolyData *output)
{
	std::cout << "ComputeLocalPCA" << std::endl;

	int NNoNorms = 0;
	int GoodNorms = 0;
	int badEvals = 0;
	int noisePoints = 0;
	int distPoints = 0;

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *Normals = vtkDoubleArray::New();
	Normals->SetNumberOfComponents(3);
	Normals->SetName("Normals");
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		double mean[3];
		double v0[3];
		double v1[3];
		double v2[3];
		double evals[3];

		if (m_Connectivity[i].size())
		{
			ComputeEigenVectorsOfPointCloud(input, m_Connectivity[i], mean, v0, v1, v2, evals);

			if ((evals[0] + evals[1] + evals[2]) > 0)
			{
				// variance explained by 3rd eigenvector
				double varexp = varexp = evals[2] / (evals[0] + evals[1] + evals[2]) * 100.0;

				if (varexp < Max3EigenValue)
				{
					double n[3];
					n[0] = v0[2];
					n[1] = v1[2];
					n[2] = v2[2];
					vtkMath::Normalize(n);

					double p[3];
					input->GetPoint(i, p);

					double pt[3];
					pt[0] = p[0] - mean[0];
					pt[1] = p[1] - mean[1];
					pt[2] = p[2] - mean[2];

					double planeDist = abs(vtkMath::Dot(pt, n));

					if (planeDist < MaxPlaneDistance)
					{
						vtkIdType id = pts->InsertNextPoint(p);
						verts->InsertNextCell(1);
						verts->InsertCellPoint(id);

						Normals->InsertNextTuple(n);

						if (ScalarMode == VTK_SCALAR_MODE_3RDEIGENVALUE)
						{
							scalars->InsertNextTuple(&varexp);
						}
						else // distance to plane
						{
							scalars->InsertNextTuple(&planeDist);
						}
						GoodNorms++;
					}
					else
					{
						distPoints++;
					}
				}
				else
				{
					noisePoints++;
				}
			}
			else
			{
				badEvals++;
			}
		}
		else
		{
			NNoNorms++;
		}
	}

	output->SetPoints(pts);
	output->SetVerts(verts);
	output->GetPointData()->SetScalars(scalars);
	output->GetPointData()->SetNormals(Normals);

	Normals->Delete();
	pts->Delete();
	verts->Delete();
	scalars->Delete();

	std::cout << "Created " << GoodNorms << " points with normals" << std::endl;
	std::cout << "Points with no neighbours to use for normal calculation " << NNoNorms << std::endl;
	std::cout << "Points with bad evals " <<  badEvals << std::endl;
	std::cout << "Points where the 3rd Eigenvector explained more than " << Max3EigenValue << " % of the variation " << noisePoints << std::endl;
	std::cout << "Points where the distance to the fitted plane was more than " << MaxPlaneDistance << " mm " << distPoints << std::endl;
}

// Generate normals for polygon meshes
int vtkPolyDataPCANormals::RequestData(
									   vtkInformation *vtkNotUsed(request),
									   vtkInformationVector **inputVector,
									   vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

//	output->DeepCopy(input);

	vtkTimerLog *timer = vtkTimerLog::New();
	timer->StartTimer();

	bool usePreConnectivity = false;
	if (usePreConnectivity)
	{
		if (SearchMode == VTK_SEARCH_MODE_SEARCHRADIUS)
		{
			CreateConnectivityUsingSearchRadius(input);
		}
		else
		{
			CreateConnectivityByPointsPerNormal(input);
		}

		ComputeLocalPCA(input, output);
	}
	else
	{
		if (SearchMode == VTK_SEARCH_MODE_SEARCHRADIUS)
		{
			ComputeLocalPCAUsingSearchRadius(input, output);
		}
		else
		{
			ComputeLocalPCAByPointsPerNormal(input, output);
		}
	}

	timer->StopTimer();
	double ElapsedTime = timer->GetElapsedTime();
	std::cout << "\nCreating PCA normals elapsedTime " << ElapsedTime << std::endl;
	timer->Delete();
	return 1;
}

void vtkPolyDataPCANormals::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

