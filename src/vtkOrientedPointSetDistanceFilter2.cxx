#include "vtkOrientedPointSetDistanceFilter2.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include <vtkPointLocator.h>
#include <iostream>
#include <vtkIntArray.h>
#include "GeneralUtils.h"
#include <vector>
#include <thread>
#include <vtkCellLocator.h>


vtkStandardNewMacro(vtkOrientedPointSetDistanceFilter2);

vtkOrientedPointSetDistanceFilter2::vtkOrientedPointSetDistanceFilter2()
{
	this->SampleSpacing = -1.0; // negative values cause the algorithm to make a reasonable guess
	DistanceMode = VTK_ORIENTEDDISTANCE_PROJECTED_MEDIAN;
	SampleFactor = 1.0;
	WeightVolume = NULL;
	WeightVolumeData = NULL;
	CreateWeightVolume  = 0;
	SearchMode = VTK_SEARCH_MODE_NNEIGHBOURS;
	NumberOfDistances = 5;
	SearchRadius = 1;
	WeightValueHigh = 3;
	ComputeMode = VTK_ORIENTEDDISTANCE_FULLVOLUME;
	InputDistances = NULL;
	MinSideVoxels = 16;
	PadVoxels = 5;
	NumThreads = 0;
}

vtkOrientedPointSetDistanceFilter2::~vtkOrientedPointSetDistanceFilter2()
{
	if (WeightVolume)
		WeightVolume->Delete();
	if (WeightVolumeData)
		WeightVolumeData->Delete();
}

int vtkOrientedPointSetDistanceFilter2::FillInputPortInformation(
  int vtkNotUsed( port ), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

int vtkOrientedPointSetDistanceFilter2::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  //// would be nice to compute the whole extent but we need more info to
  //// compute it.
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),0,1,0,1,0,1);

  // SOme of this is copied further down. 
  // This is necessary to make the class compute the output bounds correctly
  if (ComputeMode == VTK_ORIENTEDDISTANCE_FULLVOLUME || ComputeMode == VTK_ORIENTEDDISTANCE_COMBINED_FULLBAND
	  || ComputeMode == VTK_UNSIGNEDDISTANCE_FULLVOLUME)
  {
  	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
//	vtkInformation *outInfo = outputVector->GetInformationObject(0);
//	vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Get real size of point set
	double bounds[6];
	input->GetBounds(bounds);

	// The real size (probably in mm ..)
	double Xs = bounds[1] - bounds[0];
	double Ys = bounds[3] - bounds[2];
	double Zs = bounds[5] - bounds[4];

	// Find minimum sidelength
	double lmin = (bounds[1] - bounds[0]);
	int mindim = 0;
	for (int i = 0; i < 3; i++)
	{
		double l = bounds[i*2 + 1] - bounds[i*2];
		if (l < lmin)
		{
			lmin = l;
			mindim = i;
		}
	}
	int nvox = MinSideVoxels;
	SampleSpacing = lmin/nvox;

	int PadSize = PadVoxels;

	int sx = (int)(Xs/SampleSpacing + 2 * PadSize) + 1;
	int sy = (int)(Ys/SampleSpacing + 2 * PadSize) + 1;
	int sz = (int)(Zs/SampleSpacing + 2 * PadSize) + 1;

	double Xmin = bounds[0] - PadSize * SampleSpacing;
	double Ymin = bounds[2] - PadSize * SampleSpacing;
	double Zmin = bounds[4] - PadSize * SampleSpacing;

	double topleft[3] = {Xmin,Ymin,Zmin};

		// initialise the output volume
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 0, sx, 0, sy, 0, sz);
//	output->SetExtent(0, sx, 0, sy, 0, sz);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),	0, sx, 0, sy, 0, sz);
	outInfo->Set(vtkDataObject::SPACING(),	this->SampleSpacing, this->SampleSpacing, this->SampleSpacing);
	outInfo->Set(vtkDataObject::ORIGIN(), topleft, 3);
  }
  else
  {
  	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	//vtkInformation *outInfo = outputVector->GetInformationObject(0);
	//vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	double SS[3];
	InputDistances->GetSpacing(SS);
	SampleSpacing = SS[0];

	int dims[3];
	InputDistances->GetDimensions(dims);

	// Get real size of point set
	double bounds[6];
	InputDistances->GetBounds(bounds);

	double topleft[3] = {bounds[0],bounds[2],bounds[4]};

		// initialise the output volume
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
//	output->SetExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),	0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
//	output->AllocateScalars(outInfo);
	outInfo->Set(vtkDataObject::SPACING(),	this->SampleSpacing, this->SampleSpacing, this->SampleSpacing);
	outInfo->Set(vtkDataObject::ORIGIN(), topleft, 3);
  }


  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);
  return 1;
}
//
//void vtkOrientedPointSetDistanceFilter2::ExecuteInformation()
//{
//	vtkImageData *output = this->GetOutput();
//
//	output->SetScalarType(VTK_DOUBLE);
//	output->SetNumberOfScalarComponents(1);
//}

bool vtkOrientedPointSetDistanceFilter2::CalculateDistanceToSurfaceAtPosition(vtkCellLocator* locator, double* point, vtkDataSet* input, vtkDataArray* normals, vtkDoubleArray* newScalars, int offset)
{
	vtkIdType cell_id;
	int sub_id;
	double dist2, totaldist = 0;

	double tcp[3];
	locator->FindClosestPoint(point, tcp, cell_id, sub_id, dist2);

	double dist = sqrt(dist2);

	newScalars->SetValue(offset, dist);

	if (CreateWeightVolume)
	{
		double val = sqrt(vtkMath::Distance2BetweenPoints(point, tcp));

		double low = 0;
		double high = WeightValueHigh;
		double normval = 1.0 - std::max(std::min((val - low) / (high - low), 1.0), 0.0);

		WeightVolumeData->SetValue(offset, normval);
	}
	return true;
}


bool vtkOrientedPointSetDistanceFilter2::CalculateDistanceAtPosition( vtkPointLocator * locator, double * point, vtkDataSet * input, vtkDataArray * normals, vtkDoubleArray * newScalars, int offset, int LocalSearchMode )
{
	vtkIdList *neighPts = vtkIdList::New();

	if (LocalSearchMode == VTK_SEARCH_MODE_NNEIGHBOURS)
		locator->FindClosestNPoints(NumberOfDistances, point, neighPts);
	else
		locator->FindPointsWithinRadius(SearchRadius, point, neighPts);

	int NDists = neighPts->GetNumberOfIds();
	if (NDists == 0)
	{
		return false;
	}

	std::vector<double> distances;
	std::vector<int> distIDS;

	for (int n = 0; n < neighPts->GetNumberOfIds(); n++)
	{
		double cp[3];
		vtkIdType cid = neighPts->GetId(n);
		input->GetPoint(cid, cp);

		if (DistanceMode == VTK_UNSIGNEDDISTANCE_MEDIAN)
		{
			// Euclidean distance (no sign...so adopt sign from td)
			double ecd = sqrt(vtkMath::Distance2BetweenPoints(point, cp));
			double finaldist = -ecd;  // To compensate for negative value further down
		
			distances.push_back(finaldist);
			distIDS.push_back(cid);
		}
		else
		{
			// Find vector from sample point  to actual point
			double vv[3];
			vv[0] = point[0] - cp[0];
			vv[1] = point[1] - cp[1];
			vv[2] = point[2] - cp[2];

			double normal[3];
			normals->GetTuple(cid, normal);

			// projected distance (including sign)
			double td = vtkMath::Dot(vv, normal);

			double finaldist = td;
			if (DistanceMode == VTK_ORIENTEDDISTANCE_SIMPLE || DistanceMode == VTK_ORIENTEDDISTANCE_AVERAGE ||
				DistanceMode == VTK_ORIENTEDDISTANCE_MEDIAN)
			{
				// Euclidean distance (no sign...so adopt sign from td)
				double ecd = sqrt(vtkMath::Distance2BetweenPoints(point, cp));
				if (td < 0)
					finaldist = -ecd;
				else
					finaldist = ecd;
			}

			distances.push_back(finaldist);
			distIDS.push_back(cid);
		}
	}

	// Find median
	int medIndx = NDists / 2;
	CGeneralUtils::Sort2Vectors(distances, distIDS);

	// Median of distances. Id of point with median distance. Used to create weight volume
	double meddist = distances[medIndx];
	int closestID = distIDS[medIndx];

	if (DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE || DistanceMode == VTK_ORIENTEDDISTANCE_AVERAGE ||
		DistanceMode == VTK_ORIENTEDDISTANCE_SIMPLE || DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED)
	{
		double meandist = 0;
		double sdevdist = 0;
		CGeneralUtils::MeanAndSdev(distances, meandist, sdevdist);
		meddist = meandist;
	}

	newScalars->SetValue(offset, -meddist);
	neighPts->Delete();

	if (CreateWeightVolume)
	{
		double p[3];
		input->GetPoint(closestID, p);

		double val = sqrt(vtkMath::Distance2BetweenPoints(point, p));

		double low = 0;
		double high = WeightValueHigh;
		double normval = 1.0-std::max(std::min((val-low) / (high-low), 1.0), 0.0);

		WeightVolumeData->SetValue(offset, normval);
	}
	return true;
}


struct ThreadArgs
{
	int zmin;
	int zmax;
	int offsetmin;
	int offsetmax;
	int * dim;
	double * topleft;
	double SampleSpacing;
	vtkPointLocator * locator = NULL;
	vtkDataSet * input;
	vtkDataArray * normals;
	vtkDoubleArray * newScalars; 

	int SearchMode;
	int NumberOfDistances;
	double SearchRadius;
	int DistanceMode;
	int CreateWeightVolume;
	vtkDoubleArray *WeightVolumeData;
	double WeightValueHigh;
};

bool CalculateDistanceAtPositionLocalVersion( ThreadArgs Targs, double * point, int offset, int LocalSearchMode )
{
	vtkIdList *neighPts = vtkIdList::New();

	//std::cout << "Debug Findclosestpoint start. Searchmode " << LocalSearchMode << std::endl; // DEBUG
	if (LocalSearchMode == VTK_SEARCH_MODE_NNEIGHBOURS)
	{
		//std::cout << "Debug FindClosestNPoints start " << Targs.NumberOfDistances << " " << point[0] << " " << point[1] << " " << point[2] <<  std::endl; // DEBUG
		Targs.locator->FindClosestNPoints(Targs.NumberOfDistances, point, neighPts);
	}
	else
	{
		//std::cout << "Debug FindPointsWithinRadius start " << Targs.SearchRadius << " " << point[0] << " " << point[1] << " " << point[2] << std::endl; // DEBUG
		Targs.locator->FindPointsWithinRadius(Targs.SearchRadius, point, neighPts);
	}
	//std::cout << "Debug Findclosestpoint end " << std::endl;// DEBUG

	int NDists = neighPts->GetNumberOfIds();
	if (NDists == 0)
	{
		//std::cout << "Debug NDists = 0 " << std::endl;// DEBUG
		neighPts->Delete();
		return false;
	}

	std::vector<double> distances;
	std::vector<int> distIDS;

	//std::cout << "Debug ndists start " << std::endl; // DEBUG
	for (int n = 0; n < neighPts->GetNumberOfIds(); n++)
	{
		double cp[3];
		vtkIdType cid = neighPts->GetId(n);
		Targs.input->GetPoint(cid, cp);

		if (Targs.DistanceMode == VTK_UNSIGNEDDISTANCE_MEDIAN)
		{
			// Euclidean distance (no sign...so adopt sign from td)
			double ecd = sqrt(vtkMath::Distance2BetweenPoints(point, cp));
			double finaldist = -ecd;  // To compensate for negative value further down

			distances.push_back(finaldist);
			distIDS.push_back(cid);
		}
		else
		{
			// Find vector from sample point  to actual point
			double vv[3];
			vv[0] = point[0] - cp[0];
			vv[1] = point[1] - cp[1];
			vv[2] = point[2] - cp[2];

			double normal[3];
			Targs.normals->GetTuple(cid, normal);

			// projected distance (including sign)
			double td = vtkMath::Dot(vv, normal);

			double finaldist = td;
			if (Targs.DistanceMode == VTK_ORIENTEDDISTANCE_SIMPLE || Targs.DistanceMode == VTK_ORIENTEDDISTANCE_AVERAGE ||
				Targs.DistanceMode == VTK_ORIENTEDDISTANCE_MEDIAN)
			{
				// Euclidean distance (no sign...so adopt sign from td)
				double ecd = sqrt(vtkMath::Distance2BetweenPoints(point, cp));
				if (td < 0)
					finaldist = -ecd;
				else
					finaldist = ecd;
			}

			distances.push_back(finaldist);
			distIDS.push_back(cid);
		}
	}
	//std::cout << "Debug ndists end " << std::endl; // DEBUG

	// Find median
	int medIndx = NDists / 2;
	CGeneralUtils::Sort2Vectors(distances, distIDS);

	// Median of distances. Id of point with median distance. Used to create weight volume
	double meddist = distances[medIndx];
	int closestID = distIDS[medIndx];

	if (Targs.DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE || Targs.DistanceMode == VTK_ORIENTEDDISTANCE_AVERAGE ||
		Targs.DistanceMode == VTK_ORIENTEDDISTANCE_SIMPLE || Targs.DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED)
	{
		double meandist = 0;
		double sdevdist = 0;
		CGeneralUtils::MeanAndSdev(distances, meandist, sdevdist);
		meddist = meandist;
	}

	Targs.newScalars->SetValue(offset, -meddist);
	neighPts->Delete();

	if (Targs.CreateWeightVolume)
	{
		double p[3];
		Targs.input->GetPoint(closestID, p);

		double val = sqrt(vtkMath::Distance2BetweenPoints(point, p));

		double low = 0;
		double high = Targs.WeightValueHigh;
		double normval = 1.0-std::max(std::min((val-low) / (high-low), 1.0), 0.0);

		Targs.WeightVolumeData->SetValue(offset, normval);
	}
	return true;
}


void DistanceLoopThread(ThreadArgs Targs)
{
	// go through the array probing the values
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	for( z = Targs.zmin; z < Targs.zmax; z++)
	{
		std::cout << "\rCalculating distances: " << z << "/ [" << Targs.zmin << ", " << Targs.zmax << "[            ";
		zOffset = z*Targs.dim[1]*Targs.dim[0];
		point[2] = Targs.topleft[2] + z*Targs.SampleSpacing;

		for(y=0;y<Targs.dim[1];y++)
		{
			yOffset = y*Targs.dim[0] + zOffset;
			point[1] = Targs.topleft[1] + y*Targs.SampleSpacing;

			for(x=0;x<Targs.dim[0];x++)
			{
				offset = x + yOffset;
				point[0] = Targs.topleft[0] + x*Targs.SampleSpacing;

				int xtest, ytest, ztest;
				CGeneralUtils::GetXYZfromOffset(Targs.dim, offset, xtest, ytest,ztest);
				if (xtest != x || ytest != y || ztest != z)
				{
					std::cout << "something rotten with coordinate conversion" << std::endl;
				}

				int result = CalculateDistanceAtPositionLocalVersion(Targs, point, offset, Targs.SearchMode);

				// Perhaps no closest points where found
				if (!result)
				{
					result = CalculateDistanceAtPositionLocalVersion(Targs, point, offset,VTK_SEARCH_MODE_NNEIGHBOURS);
				}
			}
		}
	}
}

void DistanceLoopWithOffsetThread(ThreadArgs Targs)
{
	//std::cout << "We are in thread " << Targs.offsetmin << std::endl; // DEBUG

//	std::this_thread::sleep_for(std::chrono::milliseconds(1000000));

//		std::cout << "Starting thread " << Targs.offsetmin << std::endl; // DEBUG
	//	std::cout << "Starting thread " << Targs.offsetmin << std::endl; // DEBUG
	for (int off = Targs.offsetmin; off < Targs.offsetmax; off++)
	{
		//std::cout << "Debug print off " << off << std::endl; // DEBUG
		if (!(off % 100))
			std::cout << "\rCalculating distances: " << off << "/ [" << Targs.offsetmin << ", " << Targs.offsetmax << "[            ";

		int x,y,z;
		CGeneralUtils::GetXYZfromOffset(Targs.dim, off, x, y, z);

		double point[3];
		point[0] = Targs.topleft[0] + x*Targs.SampleSpacing;
		point[1] = Targs.topleft[1] + y*Targs.SampleSpacing;
		point[2] = Targs.topleft[2] + z*Targs.SampleSpacing;

		//std::cout << "Debug CalculateDistanceAtPositionLocalVersion before " << std::endl; // DEBUG
		int result = CalculateDistanceAtPositionLocalVersion(Targs, point, off, Targs.SearchMode);
		//std::cout << "Debug CalculateDistanceAtPositionLocalVersion after " << std::endl; // DEBUG

		// Perhaps no closest points where found
		if (!result)
		{
			//std::cout << "Debug CalculateDistanceAtPositionLocalVersion before 2" << std::endl; // DEBUG
			result = CalculateDistanceAtPositionLocalVersion(Targs, point, off,VTK_SEARCH_MODE_NNEIGHBOURS);
			//std::cout << "Debug CalculateDistanceAtPositionLocalVersion after 2" << std::endl; // DEBUG
		}
	}

	//std::cout << "We are going out of thread " << Targs.offsetmin << std::endl; // DEBUG
}

void ComputeBandedDistanceLoopWithOffsetThread(ThreadArgs Targs, vtkDoubleArray * MarkVolumeData, vtkDoubleArray * interpolatedValues)
{
	for (int off = Targs.offsetmin; off < Targs.offsetmax; off++)
	{
		//if (!(off % 100))
		//	std::cout << "\rCalculating distances: " << off << "/ [" << Targs.offsetmin << ", " << Targs.offsetmax << "[            ";

		if (MarkVolumeData->GetValue(off) == 0)
		{
			// Just copy interpolated value from previous iteration
			// and give it weight zero
			double intval = interpolatedValues->GetValue(off);
			Targs.newScalars->SetValue(off, intval);

			if (Targs.CreateWeightVolume)
			{
				Targs.WeightVolumeData->SetValue(off, 0);
			}
		}
		else
		{
			int x,y,z;
			CGeneralUtils::GetXYZfromOffset(Targs.dim, off, x, y, z);

			double point[3];
			point[0] = Targs.topleft[0] + x*Targs.SampleSpacing;
			point[1] = Targs.topleft[1] + y*Targs.SampleSpacing;
			point[2] = Targs.topleft[2] + z*Targs.SampleSpacing;


			int result = CalculateDistanceAtPositionLocalVersion(Targs, point, off, Targs.SearchMode);

			// Perhaps no closest points where found
			if (!result)
			{
				// Just copy interpolated value from previous iteration
				// and give it weight zero
				double intval = interpolatedValues->GetValue(off);
				Targs.newScalars->SetValue(off, intval);

				if (Targs.CreateWeightVolume)
				{
					Targs.WeightVolumeData->SetValue(off, 0);
				}
			}
		}
	}
}


void ComputeBandedDistanceLoopThread(ThreadArgs Targs, vtkDoubleArray * MarkVolumeData, vtkDoubleArray * interpolatedValues)
{
	// go through the array probing the values
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	for( z = Targs.zmin; z < Targs.zmax; z++)
	{
		std::cout << "\rCalculating distances: " << z << "/ [" << Targs.zmin << ", " << Targs.zmax << "[            ";
		zOffset = z*Targs.dim[1]*Targs.dim[0];
		point[2] = Targs.topleft[2] + z*Targs.SampleSpacing;

		for(y=0;y<Targs.dim[1];y++)
		{
			yOffset = y*Targs.dim[0] + zOffset;
			point[1] = Targs.topleft[1] + y*Targs.SampleSpacing;

			for(x=0;x<Targs.dim[0];x++)
			{
				offset = x + yOffset;
				point[0] = Targs.topleft[0] + x*Targs.SampleSpacing;

				if (MarkVolumeData->GetValue(offset) == 0)
				{
					// Just copy interpolated value from previous iteration
					// and give it weight zero
					double intval = interpolatedValues->GetValue(offset);
					Targs.newScalars->SetValue(offset, intval);

					if (Targs.CreateWeightVolume)
					{
						Targs.WeightVolumeData->SetValue(offset, 0);
					}
				}
				else
				{
					int result = CalculateDistanceAtPositionLocalVersion(Targs, point, offset, Targs.SearchMode);

					// Perhaps no closest points where found
					if (!result)
					{
						// Just copy interpolated value from previous iteration
						// and give it weight zero
						double intval = interpolatedValues->GetValue(offset);
						Targs.newScalars->SetValue(offset, intval);

						if (Targs.CreateWeightVolume)
						{
							Targs.WeightVolumeData->SetValue(offset, 0);
						}
					}
				}
			}
		}
	}
}



//
////void tst(int left, int right, int * dim, double * topleft, vtkPointLocator * locator) 
//void tst(int left, int right, int * dim, double * topleft, vtkDataSet * input) 
//{
//	for (int i = left; i < right; ++i) 
//	{
//		std::cout << left << right << std::endl;
//	}
//}


////Split "mem" into "parts", e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10
////if possible the function will split mem into equal chunks, if not 
////the last chunk will be slightly larger
//std::vector<int> bounds(int parts, int mem) 
//{
//	std::vector<int>bnd;
//	int delta = mem / parts;
//	int reminder = mem % parts;
//	int N1 = 0, N2 = 0;
//	bnd.push_back(N1);
//	for (int i = 0; i < parts; ++i) 
//	{
//		N2 = N1 + delta;
//		if (i == parts - 1)
//			N2 += reminder;
//		bnd.push_back(N2);
//		N1 = N2;
//	}
//	return bnd;
//}
//

//void ThreadTestHello()
//{
//	std::cout << "Hello from thread" << std::endl;
//}

void vtkOrientedPointSetDistanceFilter2::ComputeDistancesInFullVolumeMultiThreaded( vtkDataSet *input, int dim[3], vtkDoubleArray *newScalars, double topleft[3], vtkDataArray *normals )
{
	if (DistanceMode != VTK_UNSIGNEDDISTANCE_MEDIAN &&
		DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE && !normals)
	{
		std::cout << "Missing normals so no distances computed" << std::endl;
		return;
	}

	std::cout << "Computing distances in full volume (" << dim[0] * dim[1] * dim[2] << " voxels)" << std::endl;

	vtkPointLocator* locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();
	
	std::cout << "Running distance computation using " << NumThreads << " threads" << std::endl;

	const int parts = NumThreads;
//	int zmax = dim[2];
	int size = dim[0] * dim[1] * dim[2];

//	std::vector<int>bnd = CGeneralUtils::bounds(parts, zmax);
	std::vector<int>bnd = CGeneralUtils::bounds(parts, size);

	std::vector<std::thread> tt;
	////std::vector<ThreadArgs> Targs; // Make sure parm do not fall out of scope
	////Targs.resize(parts - 1);

 //   // Debug check
	//const int i = 0;
	//ThreadArgs Targs;
	//Targs.offsetmin = bnd[i];
	//Targs.offsetmax = bnd[i+1];
	//Targs.dim  = dim;
	//Targs.SampleSpacing = this->SampleSpacing;
	//Targs.topleft = topleft;
	//Targs.locator = locator;
	//Targs.input = input;
	//Targs.normals = normals;
	//Targs.newScalars = newScalars;
	//Targs.NumberOfDistances = NumberOfDistances;
	//Targs.SearchMode = SearchMode;
	//Targs.SearchRadius = SearchRadius;
	//Targs.DistanceMode = DistanceMode;
	//Targs.CreateWeightVolume = CreateWeightVolume;
	//Targs.WeightVolumeData = WeightVolumeData;
	//Targs.WeightValueHigh = WeightValueHigh;
	//std::cout << "Spawning thread " << Targs.offsetmin << std::endl; // DEBUG
	//std::thread t(DistanceLoopWithOffsetThread, Targs);
	////std::cout << "Before sleep" << std::endl;
	////std::this_thread::sleep_for(std::chrono::milliseconds(1000000));
	////std::cout << "End sleep" << std::endl;
	//t.join();


	//std::thread t(ThreadTestHello); // DEBUG Test
	//t.join();


	//std::vector<std::thread> tt;
	////std::vector<ThreadArgs> Targs; // Make sure parm do not fall out of scope
	////Targs.resize(parts - 1);

	//Launch parts-1 threads
	for (int i = 0; i < parts - 1; ++i) 
	{
		ThreadArgs Targs;
		//Targs.zmin = bnd[i];
		//Targs.zmax = bnd[i+1];
		Targs.offsetmin = bnd[i];
		Targs.offsetmax = bnd[i+1];
		Targs.dim  = dim;
		Targs.SampleSpacing = this->SampleSpacing;
		Targs.topleft = topleft;
		Targs.locator = locator;
		Targs.input = input;
		Targs.normals = normals;
		Targs.newScalars = newScalars;
		Targs.NumberOfDistances = NumberOfDistances;
		Targs.SearchMode = SearchMode;
		Targs.SearchRadius = SearchRadius;
		Targs.DistanceMode = DistanceMode;
		Targs.CreateWeightVolume = CreateWeightVolume;
		Targs.WeightVolumeData = WeightVolumeData;
		Targs.WeightValueHigh = WeightValueHigh;
		
		//std::cout << "Spawning thread " << Targs.offsetmin << std::endl; // DEBUG
 		tt.push_back(std::thread(DistanceLoopWithOffsetThread, Targs));
	}

	//Use the main thread to do part of the work !!!
	for (int i = parts - 1; i < parts; ++i) 
	{
		ThreadArgs Targs;
		//Targs.zmin = bnd[i];
		//Targs.zmax = bnd[i+1];
		Targs.offsetmin = bnd[i];
		Targs.offsetmax = bnd[i+1];
		Targs.dim  = dim;
		Targs.SampleSpacing = this->SampleSpacing;
		Targs.topleft = topleft;
		Targs.locator = locator;
		Targs.input = input;
		Targs.normals = normals;
		Targs.newScalars = newScalars;
		Targs.NumberOfDistances = NumberOfDistances;
		Targs.SearchMode = SearchMode;
		Targs.SearchRadius = SearchRadius;
		Targs.DistanceMode = DistanceMode;
		Targs.CreateWeightVolume = CreateWeightVolume;
		Targs.WeightVolumeData = WeightVolumeData;
		Targs.WeightValueHigh = WeightValueHigh;

		//std::cout << "Running extra thread in main loop " << Targs.offsetmin << std::endl; 		// DEBUG
		DistanceLoopWithOffsetThread(Targs);
	}

	//std::cout << "Before sleep" << std::endl;
	//std::this_thread::sleep_for(std::chrono::milliseconds(1000000));
	//std::cout << "End sleep" << std::endl;


	//Join parts-1 threads
	for(auto &e : tt)
	{
		//std::cout << "Joining thread" << std::endl; // DEBUG
		e.join();
	}

	if (locator)
		locator->Delete();
	if (CreateWeightVolume)
	{
		WeightVolumeData->Modified();
	}
	std::cout << std::endl;
}



void vtkOrientedPointSetDistanceFilter2::ComputeDistancesInFullVolume( vtkDataSet *input, int dim[3], vtkDoubleArray *newScalars, double topleft[3], vtkDataArray *normals )
{
	if (DistanceMode != VTK_UNSIGNEDDISTANCE_MEDIAN &&
		DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE && !normals)
	{
		std::cout << "Missing normals so no distances computed" << std::endl;
		return;
	}

	std::cout << "Computing distances in full volume (" << dim[0] * dim[1] * dim[2] << " voxels)" << std::endl;

	vtkPointLocator* locator = NULL;
	vtkCellLocator* cellLocator = NULL;
	if (DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE)
	{
		locator = vtkPointLocator::New();
		locator->SetDataSet(input);
		locator->SetNumberOfPointsPerBucket(1);
		locator->BuildLocator();
	}
	else
	{
		cellLocator = vtkCellLocator::New();
		cellLocator->SetDataSet(input);
		cellLocator->SetNumberOfCellsPerBucket(1);
		cellLocator->BuildLocator();
	}

	// go through the array probing the values
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	for(z=0;z<dim[2];z++)
	{
		std::cout << "\rCalculating distances: " << z << "/" << dim[2] << "            ";
		zOffset = z*dim[1]*dim[0];
		point[2] = topleft[2] + z*this->SampleSpacing;

		for(y=0;y<dim[1];y++)
		{
			yOffset = y*dim[0] + zOffset;
			point[1] = topleft[1] + y*this->SampleSpacing;

			for(x=0;x<dim[0];x++)
			{
				offset = x + yOffset;
				point[0] = topleft[0] + x*this->SampleSpacing;

				if (DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE)
				{
					int result = CalculateDistanceAtPosition(locator, point, input, normals, newScalars, offset, SearchMode);

					// Perhaps no closest points where found
					if (!result)
					{
						result = CalculateDistanceAtPosition(locator, point, input, normals, newScalars, offset, VTK_SEARCH_MODE_NNEIGHBOURS);
					}
				}
				else
				{
					bool result = CalculateDistanceToSurfaceAtPosition(cellLocator, point, input, normals, newScalars, offset);
				}
			}
		}
	}
	if (locator)
		locator->Delete();
	if (cellLocator)
		cellLocator->Delete();

	if (CreateWeightVolume)
	{
		WeightVolumeData->Modified();
	}
	std::cout << std::endl;
}



void vtkOrientedPointSetDistanceFilter2::BandMarkVolumeInBand( vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDoubleArray* MarkVolumeData )
{
	std::cout << "Marking narrow band with " << input->GetNumberOfPoints() << " points." << std::endl;

	int Ntot = 0;

	// Marksize
	int MS = std::max(2.0, WeightValueHigh / SampleSpacing + 1);

	std::cout << "Mark size " << MS << " voxels" << std::endl;

	// go through all cells in shape and mark narrow band
	for (int i = 0; i < input->GetNumberOfPoints(); i++)
	{
		if ((i % 500) == 0)
		{
			std::cout << "\r" << i << "    ";
		}
		double p[3];
		input->GetPoint(i, p);

		// Locate voxel where point is in
		int pbox[3];
		pbox[0] = (int)((p[0] - topleft[0]) / SampleSpacing);
		pbox[1] = (int)((p[1] - topleft[1]) / SampleSpacing);
		pbox[2] = (int)((p[2] - topleft[2]) / SampleSpacing);

		int xl = std::max(pbox[0]-MS,        0);
		int xh = std::min(pbox[0]+MS, dim[0]-1);
		int yl = std::max(pbox[1]-MS,        0);
		int yh = std::min(pbox[1]+MS, dim[1]-1);
		int zl = std::max(pbox[2]-MS,        0);
		int zh = std::min(pbox[2]+MS, dim[2]-1);


		for(int z = zl; z <= zh; z++)
		{
			int zOffset = z*dim[1]*dim[0];

			for(int y = yl; y <= yh; y++)
			{
				int yOffset = y*dim[0] + zOffset;

				for(int x = xl; x <= xh; x++)
				{
					int offset = x + yOffset;

					if (MarkVolumeData->GetValue(offset) == 0)
					{
						Ntot++;
					}
			
					MarkVolumeData->SetValue(offset, 1);
				}
			}
		}
	}
	std::cout << "\nNarrow band consists of " << Ntot << " voxels" << std::endl;
}




void vtkOrientedPointSetDistanceFilter2::ComputeDistancesInBandWithInputMultiThreaded( vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDataArray *normals )
{
	if (!InputDistances)
	{
		std::cerr << "Missing input distances" << std::endl;
		return;
	}

	std::cout << "Computing distances in band and interpolate the rest from lower resolution" << std::endl;


	// Create marking volume used for locating voxels that should be updated
	int size = dim[0] * dim[1] * dim[2];
	vtkDoubleArray * MarkVolumeData = vtkDoubleArray::New();
	MarkVolumeData->SetNumberOfComponents(1);
	MarkVolumeData->SetNumberOfValues(size);
	for (int i = 0; i < size; i++)
	{
		MarkVolumeData->SetValue(i, 0);
	}

	BandMarkVolumeInBand(input, dim, output, topleft, MarkVolumeData);

	vtkDoubleArray *newScalars = 
		vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	vtkDoubleArray *interpolatedValues = 
		vtkDoubleArray::SafeDownCast(InputDistances->GetPointData()->GetScalars());

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

	std::cout << "Running banded distance computation using " << NumThreads << " threads" << std::endl;

	const int parts = NumThreads;
	//int zmax = dim[2];
	//std::vector<int>bnd = CGeneralUtils::bounds(parts, zmax);

	std::vector<int>bnd = CGeneralUtils::bounds(parts, size);

	std::vector<std::thread> tt;

	//Launch parts-1 threads
	for (int i = 0; i < parts - 1; ++i) 
	{
		ThreadArgs Targs;
		//Targs.zmin = bnd[i];
		//Targs.zmax = bnd[i+1];
		Targs.offsetmin = bnd[i];
		Targs.offsetmax = bnd[i+1];
		Targs.dim  = dim;
		Targs.SampleSpacing = this->SampleSpacing;
		Targs.topleft = topleft;
		Targs.locator = locator;
		Targs.input = input;
		Targs.normals = normals;
		Targs.newScalars = newScalars;
		Targs.NumberOfDistances = NumberOfDistances;
		Targs.SearchMode = SearchMode;
		Targs.SearchRadius = SearchRadius;
		Targs.DistanceMode = DistanceMode;
		Targs.CreateWeightVolume = CreateWeightVolume;
		Targs.WeightVolumeData = WeightVolumeData;
		Targs.WeightValueHigh = WeightValueHigh;

		tt.push_back(std::thread(ComputeBandedDistanceLoopWithOffsetThread, Targs, MarkVolumeData, interpolatedValues));
	}

	//Use the main thread to do part of the work !!!
	for (int i = parts - 1; i < parts; ++i) 
	{
		ThreadArgs Targs;
		//Targs.zmin = bnd[i];
		//Targs.zmax = bnd[i+1];
		Targs.offsetmin = bnd[i];
		Targs.offsetmax = bnd[i+1];
		Targs.dim  = dim;
		Targs.SampleSpacing = this->SampleSpacing;
		Targs.topleft = topleft;
		Targs.locator = locator;
		Targs.input = input;
		Targs.normals = normals;
		Targs.newScalars = newScalars;
		Targs.NumberOfDistances = NumberOfDistances;
		Targs.SearchRadius = SearchRadius;
		Targs.DistanceMode = DistanceMode;
		Targs.CreateWeightVolume = CreateWeightVolume;
		Targs.WeightVolumeData = WeightVolumeData;
		Targs.WeightValueHigh = WeightValueHigh;

		ComputeBandedDistanceLoopWithOffsetThread(Targs, MarkVolumeData, interpolatedValues);
	}

	//Join parts-1 threads
	for(auto &e : tt)
	{
		e.join();
	}


//	ComputeBandedDistanceLoopThread(dim, topleft, MarkVolumeData, interpolatedValues, newScalars, locator, input, normals);

	locator->Delete();
	if (CreateWeightVolume)
	{
		WeightVolumeData->Modified();
	}
	std::cout << std::endl;
}


void vtkOrientedPointSetDistanceFilter2::ComputeDistancesInBandWithInput( vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDataArray *normals )
{
	if (!InputDistances)
	{
		std::cerr << "Missing input distances" << std::endl;
		return;
	}

	std::cout << "Computing distances in band and interpolate the rest from lower resolution" << std::endl;


	// Create marking volume used for locating voxels that should be updated
	int size = dim[0] * dim[1] * dim[2];
	vtkDoubleArray * MarkVolumeData = vtkDoubleArray::New();
	MarkVolumeData->SetNumberOfComponents(1);
	MarkVolumeData->SetNumberOfValues(size);
	for (int i = 0; i < size; i++)
	{
		MarkVolumeData->SetValue(i, 0);
	}

	BandMarkVolumeInBand(input, dim, output, topleft, MarkVolumeData);

	vtkDoubleArray *newScalars = 
		vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	vtkDoubleArray *interpolatedValues = 
		vtkDoubleArray::SafeDownCast(InputDistances->GetPointData()->GetScalars());

	vtkPointLocator* locator = NULL;
	vtkCellLocator* cellLocator = NULL;
	if (DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE)
	{
		locator = vtkPointLocator::New();
		locator->SetDataSet(input);
		locator->SetNumberOfPointsPerBucket(1);
		locator->BuildLocator();
	}
	else
	{
		cellLocator = vtkCellLocator::New();
		cellLocator->SetDataSet(input);
		cellLocator->SetNumberOfCellsPerBucket(1);
		cellLocator->BuildLocator();
	}


	//vtkPointLocator *locator = vtkPointLocator::New();
	//locator->SetDataSet(input);
	//locator->SetNumberOfPointsPerBucket(1);
	//locator->BuildLocator();

	// go through the array probing the values
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	for(z=0;z<dim[2];z++)
	{
		std::cout << "\rCalculating distances: " << z << "/" << dim[2] << "            ";
		zOffset = z*dim[1]*dim[0];
		point[2] = topleft[2] + z*this->SampleSpacing;

		for(y=0;y<dim[1];y++)
		{
			yOffset = y*dim[0] + zOffset;
			point[1] = topleft[1] + y*this->SampleSpacing;

			for(x=0;x<dim[0];x++)
			{
				offset = x + yOffset;
				point[0] = topleft[0] + x*this->SampleSpacing;

				if (MarkVolumeData->GetValue(offset) == 0)
				{
					// Just copy interpolated value from previous iteration
					// and give it weight zero
					double intval = interpolatedValues->GetValue(offset);
					newScalars->SetValue(offset, intval);

					if (CreateWeightVolume)
					{
						WeightVolumeData->SetValue(offset, 0);
					}
				}
				else
				{
					int result = false;
					if (DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE)
					{
						result = CalculateDistanceAtPosition(locator, point, input, normals, newScalars, offset, SearchMode);
					}
					else
					{
						result = CalculateDistanceToSurfaceAtPosition(cellLocator, point, input, normals, newScalars, offset);
					}

					// Perhaps no closest points where found
					if (!result)
					{
						// Just copy interpolated value from previous iteration
						// and give it weight zero
						double intval = interpolatedValues->GetValue(offset);
						newScalars->SetValue(offset, intval);

						if (CreateWeightVolume)
						{
							WeightVolumeData->SetValue(offset, 0);
						}
					}
				}
			}
		}
	}
	if (locator)
		locator->Delete();
	if (cellLocator)
		cellLocator->Delete();
	if (CreateWeightVolume)
	{
		WeightVolumeData->Modified();
	}
	std::cout << std::endl;

}


void vtkOrientedPointSetDistanceFilter2::ComputeDistancesInBand( vtkDataSet *input, int dim[3], vtkImageData* output, double topleft[3], vtkDataArray *normals )
{
	std::cout << "Computing distances in band and set the rest to 0" << std::endl;

	// Create marking volume used for locating voxels that should be updated
	int size = dim[0] * dim[1] * dim[2];
	vtkDoubleArray * MarkVolumeData = vtkDoubleArray::New();
	MarkVolumeData->SetNumberOfComponents(1);
	MarkVolumeData->SetNumberOfValues(size);
	for (int i = 0; i < size; i++)
	{
		MarkVolumeData->SetValue(i, 0);
	}

	BandMarkVolumeInBand(input, dim, output, topleft, MarkVolumeData);

	vtkDoubleArray *newScalars = 
		vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	vtkPointLocator* locator = NULL;
	vtkCellLocator* cellLocator = NULL;
	if (DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE)
	{
		locator = vtkPointLocator::New();
		locator->SetDataSet(input);
		locator->SetNumberOfPointsPerBucket(1);
		locator->BuildLocator();
	}
	else
	{
		cellLocator = vtkCellLocator::New();
		cellLocator->SetDataSet(input);
		cellLocator->SetNumberOfCellsPerBucket(1);
		cellLocator->BuildLocator();
	}


	//vtkPointLocator *locator = vtkPointLocator::New();
	//locator->SetDataSet(input);
	//locator->SetNumberOfPointsPerBucket(1);
	//locator->BuildLocator();

	// go through the array probing the values
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	for(z=0;z<dim[2];z++)
	{
		std::cout << "\rCalculating distances: " << z << "/" << dim[2] << "            ";
		zOffset = z*dim[1]*dim[0];
		point[2] = topleft[2] + z*this->SampleSpacing;

		for(y=0;y<dim[1];y++)
		{
			yOffset = y*dim[0] + zOffset;
			point[1] = topleft[1] + y*this->SampleSpacing;

			for(x=0;x<dim[0];x++)
			{
				offset = x + yOffset;
				point[0] = topleft[0] + x*this->SampleSpacing;

				if (MarkVolumeData->GetValue(offset) == 0)
				{
					// set value to 0 and give it weight zero
					newScalars->SetValue(offset, 0);

					if (CreateWeightVolume)
					{
						WeightVolumeData->SetValue(offset, 0);
					}
				}
				else
				{
					int result = false;
					if (DistanceMode != VTK_UNSIGNEDDISTANCE_SURFACE)
					{
						result = CalculateDistanceAtPosition(locator, point, input, normals, newScalars, offset, SearchMode);
					}
					else
					{
						result = CalculateDistanceToSurfaceAtPosition(cellLocator, point, input, normals, newScalars, offset);
					}

					// Perhaps no closest points where found
					if (!result)
					{
						// set value to 0 and give it weight zero
						newScalars->SetValue(offset, 0);

						if (CreateWeightVolume)
						{
							WeightVolumeData->SetValue(offset, 0);
						}
					}
				}
			}
		}
	}
	if (locator)
		locator->Delete();
	if (cellLocator)
		cellLocator->Delete();
	if (CreateWeightVolume)
	{
		WeightVolumeData->Modified();
	}
	std::cout << std::endl;
}

//void vtkOrientedPointSetDistanceFilter2::CreateNewVolumeFromPolyData()
//{
//	vtkDataSet *input = this->GetInput();
//
//	// Get real size of point set
//	double bounds[6];
//	input->GetBounds(bounds);
//
//	// The real size (probably in mm ..)
//	double Xs = bounds[1] - bounds[0];
//	double Ys = bounds[3] - bounds[2];
//	double Zs = bounds[5] - bounds[4];
//
//	// Find minimum sidelength
//	double lmin = (bounds[1] - bounds[0]);
//	int mindim = 0;
//	for (int i = 0; i < 3; i++)
//	{
//		double l = bounds[i*2 + 1] - bounds[i*2];
//		if (l < lmin)
//		{
//			lmin = l;
//			mindim = i;
//		}
//	}
//	int nvox = MinSideVoxels;
//	SampleSpacing = lmin/nvox;
//
//	int PadSize = PadVoxels;
//
//	int sx = (int)(Xs/SampleSpacing + 2 * PadSize) + 1;
//	int sy = (int)(Ys/SampleSpacing + 2 * PadSize) + 1;
//	int sz = (int)(Zs/SampleSpacing + 2 * PadSize) + 1;
//
//	double Xmin = bounds[0] - PadSize * SampleSpacing;
//	double Ymin = bounds[2] - PadSize * SampleSpacing;
//	double Zmin = bounds[4] - PadSize * SampleSpacing;
//
//	double topleft[3] = {Xmin,Ymin,Zmin};
//
//	this->GetOutput()->SetWholeExtent(0, sx, 0, sy, 0, sz);
//	this->GetOutput()->SetUpdateExtent(0, sx, 0, sy, 0, sz);
//	this->GetOutput()->SetOrigin(topleft);
//	this->GetOutput()->SetSpacing(SampleSpacing, SampleSpacing, SampleSpacing);
//}

void vtkOrientedPointSetDistanceFilter2::CreateNewVolumeFromPolyData(vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Get real size of point set
	double bounds[6];
	input->GetBounds(bounds);

	// The real size (probably in mm ..)
	double Xs = bounds[1] - bounds[0];
	double Ys = bounds[3] - bounds[2];
	double Zs = bounds[5] - bounds[4];

	// Find minimum sidelength
	double lmin = (bounds[1] - bounds[0]);
	int mindim = 0;
	for (int i = 0; i < 3; i++)
	{
		double l = bounds[i*2 + 1] - bounds[i*2];
		if (l < lmin)
		{
			lmin = l;
			mindim = i;
		}
	}
	int nvox = MinSideVoxels;
	SampleSpacing = lmin/nvox;

	int PadSize = PadVoxels;

	int sx = (int)(Xs/SampleSpacing + 2 * PadSize) + 1;
	int sy = (int)(Ys/SampleSpacing + 2 * PadSize) + 1;
	int sz = (int)(Zs/SampleSpacing + 2 * PadSize) + 1;

	double Xmin = bounds[0] - PadSize * SampleSpacing;
	double Ymin = bounds[2] - PadSize * SampleSpacing;
	double Zmin = bounds[4] - PadSize * SampleSpacing;

	double topleft[3] = {Xmin,Ymin,Zmin};

		// initialise the output volume
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 0, sx, 0, sy, 0, sz);
	output->SetExtent(0, sx, 0, sy, 0, sz);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),	0, sx, 0, sy, 0, sz);
	output->AllocateScalars(outInfo);
	outInfo->Set(vtkDataObject::SPACING(),	this->SampleSpacing, this->SampleSpacing, this->SampleSpacing);
	outInfo->Set(vtkDataObject::ORIGIN(), topleft, 3);
}

void vtkOrientedPointSetDistanceFilter2::CreateNewVolumeFromInputVolume(vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	if (InputDistances == NULL)
	{
		std::cerr << "Problem with input distances" << std::endl;
		return;
	}

	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	double SS[3];
	InputDistances->GetSpacing(SS);
	SampleSpacing = SS[0];

	int dims[3];
	InputDistances->GetDimensions(dims);

	// Get real size of point set
	double bounds[6];
	InputDistances->GetBounds(bounds);

	double topleft[3] = {bounds[0],bounds[2],bounds[4]};

		// initialise the output volume
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
	output->SetExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),	0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
	output->AllocateScalars(outInfo);
	outInfo->Set(vtkDataObject::SPACING(),	this->SampleSpacing, this->SampleSpacing, this->SampleSpacing);
	outInfo->Set(vtkDataObject::ORIGIN(), topleft, 3);
}

int vtkOrientedPointSetDistanceFilter2::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the input
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // get the output
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


	if (ComputeMode == VTK_ORIENTEDDISTANCE_FULLVOLUME || ComputeMode == VTK_ORIENTEDDISTANCE_COMBINED_FULLBAND ||
		ComputeMode == VTK_UNSIGNEDDISTANCE_FULLVOLUME)
	{
		CreateNewVolumeFromPolyData(inputVector, outputVector);
	}
	else
	{
		CreateNewVolumeFromInputVolume(inputVector, outputVector);
	}

	vtkDoubleArray *newScalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	double bounds[6];
	output->GetBounds(bounds);
	int dim[3];
	output->GetDimensions(dim);

	double topleft[3] = {bounds[0],bounds[2],bounds[4]};

	std::cout <<"Created output volume of dimensions: ("
		<< dim[0] << ", " << dim[1] << ", " << dim[2] << ") size: " << dim[0] * dim[1] * dim[2]  << std::endl;
	std::cout << "Samplespacing " << SampleSpacing << std::endl;

	vtkDebugMacro(<<"Created output volume of dimensions: ("
		<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" );

	vtkDataArray *normals = input->GetPointData()->GetNormals(); 

	if (!normals)
	{
		std::cerr << "No normals found" << std::endl;
	}
	//else
	//{
	//	std::cout << "Number of normals " << normals->GetNumberOfTuples() << std::endl;
	//}

	if (CreateWeightVolume)
	{
		// Update weight value so it fits the resolution
		WeightValueHigh = std::max(WeightValueHigh, SampleSpacing * 2);

		std::cout << "Using weight value high: " << WeightValueHigh << std::endl;

		int size = dim[0] * dim[1] * dim[2];
		WeightVolumeData = vtkDoubleArray::New();
		WeightVolumeData->SetNumberOfComponents(1);
		WeightVolumeData->SetNumberOfValues(size);

		WeightVolume = vtkImageData::New();
		WeightVolume->CopyStructure(output);
		WeightVolume->GetPointData()->SetScalars(WeightVolumeData);
	}

	std::cout << input->GetNumberOfPoints() << " points in input" << std::endl;

	if (DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE)
		std::cout << "Computing distances using average of projected distances" << std::endl;
	else if (DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED_MEDIAN)
		std::cout << "Computing distances using median of projected distances" << std::endl;
	else if (DistanceMode == VTK_ORIENTEDDISTANCE_MEDIAN)
		std::cout << "Computing distances using median of Euclidean distances" << std::endl;
	else if (DistanceMode == VTK_ORIENTEDDISTANCE_AVERAGE)
		std::cout << "Computing distances using average of Euclidean distances" << std::endl;
	else if (DistanceMode == VTK_UNSIGNEDDISTANCE_MEDIAN)
		std::cout << "Computing unsigned distances using median of Euclidean distances" << std::endl;
	else if (DistanceMode == VTK_UNSIGNEDDISTANCE_SURFACE)
		std::cout << "Computing unsigned distances to surface" << std::endl;

	if (SearchMode == VTK_SEARCH_MODE_NNEIGHBOURS)
		std::cout << "Using nearest neighbour search mode with " << NumberOfDistances << " neighbours to be found" << std::endl;
	else if (SearchMode == VTK_SEARCH_MODE_SEARCHRADIUS)
		std::cout << "Using search radius of " << SearchRadius << " to locate neighbours" << std::endl;
	else 
		std::cout << "Unknown search mode chosen?" << std::endl;

	// NumThreads = 1; 
	if (NumThreads != -1)
	{
		// If NumThreads = 0 we automatically determine the number of cores. We leave one core free
		if (NumThreads == 0)
		{	
			NumThreads = std::thread::hardware_concurrency() - 1;
		}
	}
	if (DistanceMode == VTK_UNSIGNEDDISTANCE_SURFACE)
	{
		// vtkCellLocator is not thread safe
		NumThreads = 1;
	}

	if (ComputeMode == VTK_ORIENTEDDISTANCE_FULLVOLUME || ComputeMode == VTK_UNSIGNEDDISTANCE_FULLVOLUME)
	{
		// If less than two threads we just use single mode
		if (NumThreads <= 1)
			ComputeDistancesInFullVolume(input, dim, newScalars, topleft, normals);	
		else
			ComputeDistancesInFullVolumeMultiThreaded(input, dim, newScalars, topleft, normals);
	}
	else if (ComputeMode == VTK_ORIENTEDDISTANCE_BAND)
	{
		// If less than two threads we just use single mode
		if (NumThreads <= 1)
			ComputeDistancesInBandWithInput(input, dim, output, topleft, normals);	
		else
			ComputeDistancesInBandWithInputMultiThreaded(input, dim, output, topleft, normals);	
	}
	else if (ComputeMode == VTK_ORIENTEDDISTANCE_COMBINED_FULLBAND)
	{
		ComputeDistancesInBand(input, dim, output, topleft, normals);	
	}
	return 1;
}

void vtkOrientedPointSetDistanceFilter2::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);

	os << indent << "Sample Spacing:" << this->SampleSpacing << "\n";
}

