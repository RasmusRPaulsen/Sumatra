#include "vtkPolyDataOrientNormalsByVoting.h"

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

#include <vtkExtMisc.h>
#include <vtkDoubleArray.h>
#include <vtkPointLocator.h>

#include <deque>

vtkStandardNewMacro(vtkPolyDataOrientNormalsByVoting);


vtkPolyDataOrientNormalsByVoting::vtkPolyDataOrientNormalsByVoting()
{
	SearchRadius = 1.0;
	MinCliqueSize = 5;
	NeighbourAngle = 15;
	LargestCliqueOnly = false;
	ReferencePolyData = NULL;
}

void vtkPolyDataOrientNormalsByVoting::CreateConnectivity(vtkPolyData *input)
{
	std::cout << "Creating connectivity" << std::endl;

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


// Generate normals for polygon meshes
int vtkPolyDataOrientNormalsByVoting::RequestData(
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

	if (input->GetNumberOfPoints() == 0)
	{
		std::cerr << "No points in input" << std::endl;
		return 0;
	}

	output->DeepCopy(input);

	vtkTimerLog *timer = vtkTimerLog::New();
	timer->StartTimer();

	bool UseConnectivityGraph = false;
	if (UseConnectivityGraph)
	{
		CreateConnectivity(input);
		ComputeCliques(input, output);
	}
	else
	{
		ComputeCliquesNoConnectivity(input, output);
	}

	timer->StopTimer();
	double ElapsedTime = timer->GetElapsedTime();
	std::cout << "\nOrient normals by voting took " << ElapsedTime << std::endl;
	timer->Delete();

	return 1;
}

bool vtkPolyDataOrientNormalsByVoting::PointsConnected(int pid1, int pid2, vtkPolyData *input, vtkDataArray *normals, bool &flip)
{
	double MaxAngle = NeighbourAngle;
//	double MaxAngleR = cos(30 * vtkMath::DegreesToRadians());

	double n1[3];
	double n2[3];

	normals->GetTuple(pid1, n1);
	normals->GetTuple(pid2, n2);

//	double angle = acos(vtkMath::Dot(n1, n2)) * vtkMath::RadiansToDegrees();
	double angle = vtkMath::DegreesFromRadians(acos(vtkMath::Dot(n1, n2)));

	//if (abs(vtkMath::Dot(n1, n2)) > MaxAngleR)
	//{
	//	return false;
	//}

	flip = false;

	if (abs(angle) > 90)
	{
		flip = true;
		angle = abs(angle-180);
	}

	if (abs(angle) > MaxAngle)
	{
		return false;
	}
	
	return true;
}

void vtkPolyDataOrientNormalsByVoting::KeepLargestClique( vtkPolyData *input, vtkPolyData *output )
{
	// Simple copy points from the larger cliques
	// by creating a new polydata
	vtkDataArray *normals = output->GetPointData()->GetNormals();
	if (!normals)
	{
		std::cerr << "No normals" << std::endl;
		return;
	}

	vtkPolyData *pd = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *outnormals = vtkDoubleArray::New();
	outnormals->SetNumberOfComponents(3);
	outnormals->SetName("Normals");

	// Find largest clique
	size_t maxSize = 0;
	int maxID   = -1;
	for (unsigned int i = 0; i < m_Cliques.size(); i++)
	{
		if (m_Cliques[i].size() > maxSize)
		{
			maxSize = m_Cliques[i].size();
			maxID = i;
		}
	}

	if (maxID == -1)
	{
		std::cerr << "No max clique found?" << std::endl;
		return;
	}

	for (unsigned int j = 0; j < m_Cliques[maxID].size(); j++)
	{
		int pid = m_Cliques[maxID][j];
		double p[3];
		double n[3];

		output->GetPoint(pid, p);
		normals->GetTuple(pid, n);

		vtkIdType id = pts->InsertNextPoint(p);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);

		outnormals->InsertNextTuple(n);
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->SetNormals(outnormals);

	outnormals->Delete();
	pts->Delete();
	verts->Delete();

	//Finally copy result
	output->DeepCopy(pd);

	pd->Delete();

	std::cout << "Kept largest clique with " << m_Cliques[maxID].size() << " points" << std::endl;
}

void vtkPolyDataOrientNormalsByVoting::CutCliques( vtkPolyData *input, vtkPolyData *output )
{
	// Simple copy points from the larger cliques
	// by creating a new polydata

	vtkDataArray *normals = output->GetPointData()->GetNormals();
	if (!normals)
	{
		std::cerr << "No normals" << std::endl;
		return;
	}
	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
	if (!scalars)
	{
		std::cerr << "No scalars" << std::endl;
		return;
	}

	vtkPolyData *pd = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *outnormals = vtkDoubleArray::New();
	outnormals->SetNumberOfComponents(3);
	outnormals->SetName("Normals");
	vtkDoubleArray *outscalars = vtkDoubleArray::New();
	outscalars->SetNumberOfComponents(1);

	int cliqCopy = 0;
	int cliqCut  = 0;
	for (unsigned int i = 0; i < m_Cliques.size(); i++)
	{
		if ((int)m_Cliques[i].size() > MinCliqueSize)
		{
			for (unsigned int j = 0; j < m_Cliques[i].size(); j++)
			{
				int pid = m_Cliques[i][j];
				double p[3];
				double n[3];
				double s = cliqCopy;

				output->GetPoint(pid, p);
				normals->GetTuple(pid, n);
//				s = scalars->GetValue(pid);

				vtkIdType id = pts->InsertNextPoint(p);
				verts->InsertNextCell(1);
				verts->InsertCellPoint(id);

				outnormals->InsertNextTuple(n);
				outscalars->InsertNextTuple(&s);
			}
			cliqCopy++;
		}
		else
		{
			cliqCut++;
		}
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->SetScalars(outscalars);
	pd->GetPointData()->SetNormals(outnormals);

	outnormals->Delete();
	pts->Delete();
	verts->Delete();
	outscalars->Delete();

	//Finally copy result
	output->DeepCopy(pd);

	pd->Delete();

	std::cout << "Cliques copied " << cliqCopy << " and cut " << cliqCut << std::endl;
}


void vtkPolyDataOrientNormalsByVoting::ComputeCliques( vtkPolyData *input, vtkPolyData *output )
{
	std::cout << "Computing cliques" << std::endl;

	vtkDataArray *normals = input->GetPointData()->GetNormals();
	if (!normals)
	{
		std::cerr << "No normals" << std::endl;
		return;
	}

	int numPoints = input->GetNumberOfPoints();

	// 0 = not visited
	// 1 = visited and part of clique
	std::vector<int> status(numPoints, 0);
	
	int CliqueNum = 1;

	for (int i = 0; i < numPoints; i++)
	{
		if (status[i] == 0)
		{
			std::deque<vtkIdType> visitlist;

			visitlist.push_back(i);

			std::vector<int> clique;

			clique.push_back(i);

			while (!visitlist.empty())
			{
				vtkIdType curId = visitlist.front(); 
				visitlist.pop_front();

				status[curId] = CliqueNum;

				for (unsigned int nc = 0; nc < m_Connectivity[curId].size(); nc++)
				{
					vtkIdType cid = m_Connectivity[curId][nc];

					if (status[cid] == 0)
					{
						bool flip = false;
						if (PointsConnected(curId, cid, input, normals, flip))
						{
							status[cid] = CliqueNum;
							visitlist.push_back(cid);
							clique.push_back(cid);

							if (flip)
							{
								double nn[3];
								normals->GetTuple(cid, nn);
								nn[0] = -nn[0];
								nn[1] = -nn[1];
								nn[2] = -nn[2];
								normals->SetTuple(cid, nn);
							}
						}
					}
				}
			}
			if (clique.size() > 1)
			{
				m_Cliques.push_back(clique);
				CliqueNum++;
			}
		}
	}
	std::cout << "Found " << m_Cliques.size() << " cliques" << std::endl;

	// Hack to make it accept changes
	normals->Modified();
	output->GetPointData()->SetNormals(normals);

	// Assign clique numbers to point scalars
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < output->GetNumberOfPoints(); i++)
	{
		double val = (double)status[i];
		scalars->InsertNextTuple(&val);
	}
	output->GetPointData()->SetScalars(scalars);
	scalars->Delete();

	if (ReferencePolyData != NULL)
	{
		if (!UseReferencePolyDataInVoting(input, output))
		{
			MaxSpanningInVoting(input, output);
		}
	}
	else
	{
		MaxSpanningInVoting(input, output);
	}

	if (MinCliqueSize > 0)
	{
		if (LargestCliqueOnly)
		{
			KeepLargestClique(input, output);
		}
		else
		{
			CutCliques(input, output);
		}
	}
}

void vtkPolyDataOrientNormalsByVoting::ComputeCliquesNoConnectivity( vtkPolyData *input, vtkPolyData *output )
{
	std::cout << "Computing cliques using a search radius of " << SearchRadius << std::endl;

	vtkDataArray *normals = input->GetPointData()->GetNormals();
	if (!normals)
	{
		std::cerr << "No normals" << std::endl;
		return;
	}

	vtkPointLocator *Locator = vtkPointLocator::New();
	Locator->SetDataSet(input);
	Locator->SetNumberOfPointsPerBucket(1);
	Locator->BuildLocator();

	int numPoints = input->GetNumberOfPoints();

	// 0 = not visited
	// 1 = visited and part of clique
	std::vector<int> status(numPoints, 0);

	// Starting clique number
	int CliqueNum = 1;

	// Compute some statistics
	int NNeighbourLookups = 0;
	int TotalNeighbourfound = 0;

	for (int i = 0; i < numPoints; i++)
	{
		if (status[i] == 0)
		{
			std::deque<vtkIdType> visitlist;
			visitlist.push_back(i);

			// New clique
			std::vector<int> clique;

			clique.push_back(i);

			while (!visitlist.empty())
			{
				vtkIdType curId = visitlist.front(); 
				visitlist.pop_front();

				status[curId] = CliqueNum;

				double p[3];
				input->GetPoint(curId, p);

				vtkIdList *idlist = vtkIdList::New();
				Locator->FindPointsWithinRadius(SearchRadius, p, idlist);

				NNeighbourLookups++;
				TotalNeighbourfound += idlist->GetNumberOfIds();

				for (int nc = 0; nc < idlist->GetNumberOfIds(); nc++)
				{
					vtkIdType cid = idlist->GetId(nc);

					if (status[cid] == 0)
					{
						bool flip = false;
						if (PointsConnected(curId, cid, input, normals, flip))
						{
							status[cid] = CliqueNum;
							visitlist.push_back(cid);
							clique.push_back(cid);

							if (flip)
							{
								double nn[3];
								normals->GetTuple(cid, nn);
								nn[0] = -nn[0];
								nn[1] = -nn[1];
								nn[2] = -nn[2];
								normals->SetTuple(cid, nn);
							}
						}
					}
				}
				idlist->Delete();
			}
			if (clique.size() > 1)
			{
				m_Cliques.push_back(clique);
				CliqueNum++;
			}
		}
	}
	std::cout << "Found " << m_Cliques.size() << " cliques. Average neighbours per neighbour lookup " << (double)TotalNeighbourfound / (double)NNeighbourLookups 
		<< " lookups: " << NNeighbourLookups << std::endl;

	// Hack to make it accept changes
	normals->Modified();
	output->GetPointData()->SetNormals(normals);

	// Assign clique numbers to point scalars
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < output->GetNumberOfPoints(); i++)
	{
		double val = (double)status[i];
		scalars->InsertNextTuple(&val);
	}
	output->GetPointData()->SetScalars(scalars);
	scalars->Delete();

	if (ReferencePolyData != NULL)
	{
		if (!UseReferencePolyDataInVoting(input, output))
		{
			MaxSpanningInVoting(input, output);
		}
	}
	else
	{
		MaxSpanningInVoting(input, output);
	}

	if (MinCliqueSize > 0)
	{
		if (LargestCliqueOnly)
		{
			KeepLargestClique(input, output);
		}
		else
		{
			CutCliques(input, output);
		}
	}
	Locator->Delete();
}


bool vtkPolyDataOrientNormalsByVoting::UseReferencePolyDataInVoting( vtkPolyData *input, vtkPolyData *output )
{
	vtkDataArray *normals = input->GetPointData()->GetNormals();
	if (!normals)
	{
		std::cerr << "No normals" << std::endl;
		return false;
	}

	vtkDataArray *refnorms = ReferencePolyData->GetPointData()->GetNormals();
	if (!refnorms)
	{
		std::cerr << "No normals defined for reference poly data" << std::endl;
		return false;
	}

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(ReferencePolyData);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();


	for (unsigned int i = 0; i < m_Cliques.size(); i++)
	{
		int goodnorms = 0;
		int badnorms = 0;
		for (unsigned int j = 0; j < m_Cliques[i].size(); j++)
		{
			// Reference id
			int pid = m_Cliques[i][j];
			double p[3];
			double n[3];

			output->GetPoint(pid, p);
			normals->GetTuple(pid, n);

			vtkIdType cid = locator->FindClosestPoint(p);
			double refn[3];
			double refp[3];
			ReferencePolyData->GetPoint(cid, refp);
			refnorms->GetTuple(cid, refn);

			double dist2 = vtkMath::Distance2BetweenPoints(refp, p);

			bool badn = (vtkMath::Dot(refn, n) < 0);
			if (badn)
			{
				badnorms++;
			}
			else
			{
				goodnorms++;
			}
		}
		bool flip = (goodnorms < badnorms);
//		std::cout << "Clique " << i << " normals: good " << goodnorms << " bad " << badnorms << " flip: " << flip << std::endl;

		if (flip)
		{
			for (unsigned int j = 0; j < m_Cliques[i].size(); j++)
			{
				// Reference id
				int pid = m_Cliques[i][j];
				double n[3];

				normals->GetTuple(pid, n);

				n[0] = -n[0];
				n[1] = -n[1];
				n[2] = -n[2];

				normals->SetTuple(pid, n);
			}
		}
	}
	normals->Modified();
	output->GetPointData()->SetNormals(normals);
	locator->Delete();

	return true;
}

void vtkPolyDataOrientNormalsByVoting::MaxSpanningInVoting( vtkPolyData *input, vtkPolyData *output )
{
	vtkDataArray *normals = input->GetPointData()->GetNormals();
	if (!normals)
	{
		std::cerr << "No normals" << std::endl;
		return;
	}

	for (unsigned int i = 0; i < m_Cliques.size(); i++)
	{

		vtkPoints *refPoints = vtkPoints::New();

		// Normal end points - should always point outwards
		// Therefore the mean distance to the center of mass should be bigger for them
		vtkPoints *newPoints = vtkPoints::New();

		for (unsigned int j = 0; j < m_Cliques[i].size(); j++)
		{
			// Reference id
			int pid = m_Cliques[i][j];
			double p[3];
			double n[3];

			output->GetPoint(pid, p);
			normals->GetTuple(pid, n);

			double np[3];
			np[0] = p[0] + n[0];
			np[1] = p[1] + n[1];
			np[2] = p[2] + n[2];

			newPoints->InsertNextPoint(np);
			refPoints->InsertNextPoint(p);
		}
		
		double DistMeanRef = 0;
		double DistMeanNormend = 0;
		double DistSdev = 0;
		vtkExtMisc::DistanceFromCMStats(refPoints, DistMeanRef, DistSdev);
		vtkExtMisc::DistanceFromCMStats(newPoints, DistMeanNormend, DistSdev);

		newPoints->Delete();
		refPoints->Delete();

		bool flip = (DistMeanNormend > DistMeanRef);

		if (flip)
		{
			for (unsigned int j = 0; j < m_Cliques[i].size(); j++)
			{
				// Reference id
				int pid = m_Cliques[i][j];
				double n[3];

				normals->GetTuple(pid, n);

				n[0] = -n[0];
				n[1] = -n[1];
				n[2] = -n[2];

				normals->SetTuple(pid, n);
			}
		}
	}
	normals->Modified();
	output->GetPointData()->SetNormals(normals);
}


void vtkPolyDataOrientNormalsByVoting::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

