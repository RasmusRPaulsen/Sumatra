#include "vtkOrientedPointSetDistanceFilter.h"

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

vtkStandardNewMacro(vtkOrientedPointSetDistanceFilter);

vtkOrientedPointSetDistanceFilter::vtkOrientedPointSetDistanceFilter()
{
	this->SampleSpacing = -1.0; // negative values cause the algorithm to make a reasonable guess
	DistanceMode = VTK_ORIENTEDDISTANCE_SIMPLE;
	SampleFactor = 1.0;
	ReferenceVolume = NULL;
	ReferenceVolumeData = NULL;
	CreateReferenceVolume  = 0;
	NumberOfDistances = 5;
}

vtkOrientedPointSetDistanceFilter::~vtkOrientedPointSetDistanceFilter()
{
	if (ReferenceVolume)
		ReferenceVolume->Delete();
	if (ReferenceVolumeData)
		ReferenceVolumeData->Delete();
}

int vtkOrientedPointSetDistanceFilter::FillInputPortInformation(
  int vtkNotUsed( port ), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

int vtkOrientedPointSetDistanceFilter::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // would be nice to compute the whole extent but we need more info to
  // compute it.
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),0,1,0,1,0,1);

  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);
  return 1;
}
//
//
//void vtkOrientedPointSetDistanceFilter::ExecuteInformation()
//{
//	vtkImageData *output = this->GetOutput();
//
//	output->SetScalarType(VTK_DOUBLE);
//	output->SetNumberOfScalarComponents(1);
//
//	// would be nice to compute the whole extent but we need more info to
//	// compute it.
//	//	output->SetWholeExtent(0,0,0,0,0,0);
//}


void vtkOrientedPointSetDistanceFilter::UseSimpleCloseNeighbour(vtkDataSet *input, int dim[3],
																vtkDoubleArray *newScalars, double topleft[3],
																vtkDataArray *normals)
{
	std::cout << "Computing distances using Euclidean distances " << std::endl;

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

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

				// closest point
				double cp[3];
				vtkIdType cid = locator->FindClosestPoint(point);

				input->GetPoint(cid, cp);

				double dist2 = vtkMath::Distance2BetweenPoints(point, cp);

				double td = sqrt(dist2);

				// Find vector from point on surface to actual point
				double vv[3];
				vv[0] = point[0] - cp[0];
				vv[1] = point[1] - cp[1];
				vv[2] = point[2] - cp[2];

				double factor = 1;

				if (normals)
				{
					double normal[3];

					normals->GetTuple(cid, normal);

					factor = -vtkMath::Dot(vv, normal);

					// Normalise
					factor /= fabs(factor);
				}

				newScalars->SetValue(offset, factor*td);
			}
		}
	}
	locator->Delete();
	std::cout << std::endl;
}

void vtkOrientedPointSetDistanceFilter::UseProjectedCloseNeighbour(vtkDataSet *input, int dim[3],
																   vtkDoubleArray *newScalars, double topleft[3],
																   vtkDataArray *normals)
{
	if (!normals)
		return;

	std::cout << "Computing distances using projected distances " << std::endl;

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

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

				// closest point
				double cp[3];
				vtkIdType cid = locator->FindClosestPoint(point);

				input->GetPoint(cid, cp);

				// Find vector from point on surface to actual point
				double vv[3];
				vv[0] = point[0] - cp[0];
				vv[1] = point[1] - cp[1];
				vv[2] = point[2] - cp[2];

				double normal[3];
				normals->GetTuple(cid, normal);

				double td = vtkMath::Dot(vv, normal);
				newScalars->SetValue(offset, -td);
				if (CreateReferenceVolume)
				{
					ReferenceVolumeData->SetValue(offset, cid);
				}

				//vtkIdList *neighPts = vtkIdList::New();
				//locator->FindClosestNPoints(1, point, neighPts);
				//
				//double mindist = VTK_LARGE_FLOAT; // dummy value
				//int closestID = -1;

				//for (int n = 0; n < neighPts->GetNumberOfIds(); n++)
				//{
				//	vtkIdType cid = neighPts->GetId(n);

				//	input->GetPoint(cid, cp);

				//	double td = 0;

				//	// Find vector from point on surface to actual point
				//	double vv[3];
				//	vv[0] = point[0] - cp[0];
				//	vv[1] = point[1] - cp[1];
				//	vv[2] = point[2] - cp[2];

				//	double normal[3];
				//	normals->GetTuple(cid, normal);

				//	td = vtkMath::Dot(vv, normal);
				//	if (abs(td) < abs(mindist))
				//	{
				//		mindist = td;
				//		closestID = cid;
				//	}
				//}
				//newScalars->SetValue(offset, -mindist);
				//neighPts->Delete();
				//if (CreateReferenceVolume)
				//{
				//	ReferenceVolume->SetValue(offset, closestID);
				//}
			}
		}
	}
	locator->Delete();
	if (CreateReferenceVolume)
	{
		ReferenceVolumeData->Modified();
	}
}


void vtkOrientedPointSetDistanceFilter::UseProjectedCloseNeighbourMedian( vtkDataSet *input, int dim[3], vtkDoubleArray *newScalars, double topleft[3], vtkDataArray *normals )
{
	if (!normals)
		return;

	std::cout << "Computing distances using median of projected distances " << std::endl;

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

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

				vtkIdList *neighPts = vtkIdList::New();
				locator->FindClosestNPoints(NumberOfDistances, point, neighPts);
				
				// Keep track of both distances and their ID, so they can be used in a lookup later
				std::vector<double> distances;
				std::vector<int> distIDS;

				for (int n = 0; n < neighPts->GetNumberOfIds(); n++)
				{
					double cp[3];

					vtkIdType cid = neighPts->GetId(n);

					input->GetPoint(cid, cp);

					double td = 0;

					// Find vector from point on surface to actual point
					double vv[3];
					vv[0] = point[0] - cp[0];
					vv[1] = point[1] - cp[1];
					vv[2] = point[2] - cp[2];

					double normal[3];
					normals->GetTuple(cid, normal);

					td = vtkMath::Dot(vv, normal);
					distances.push_back(td);
					distIDS.push_back(cid);
				}
				int medIndx = NumberOfDistances / 2;
				CGeneralUtils::Sort2Vectors(distances, distIDS);

				double meddist = distances[medIndx];
				int closestID = distIDS[medIndx];

				newScalars->SetValue(offset, -meddist);
				neighPts->Delete();
				if (CreateReferenceVolume)
				{
					ReferenceVolumeData->SetValue(offset, closestID);
				}
			}
		}
	}
	locator->Delete();
	if (CreateReferenceVolume)
	{
		ReferenceVolumeData->Modified();
	}
	std::cout << std::endl;
}


void vtkOrientedPointSetDistanceFilter::UseProjectedCloseNeighbourAverage( vtkDataSet *input, int dim[3], vtkDoubleArray *newScalars, double topleft[3], vtkDataArray *normals )
{
	if (!normals)
		return;

	std::cout << "Computing distances using average of projected distances " << std::endl;

	vtkPointLocator *locator = vtkPointLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfPointsPerBucket(1);
	locator->BuildLocator();

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

				vtkIdList *neighPts = vtkIdList::New();
				locator->FindClosestNPoints(NumberOfDistances, point, neighPts);

				// Keep track of both distances and their ID, so they can be used in a lookup later
				std::vector<double> distances;
				std::vector<int> distIDS;

				for (int n = 0; n < neighPts->GetNumberOfIds(); n++)
				{
					double cp[3];

					vtkIdType cid = neighPts->GetId(n);

					input->GetPoint(cid, cp);

					double td = 0;

					// Find vector from point on surface to actual point
					double vv[3];
					vv[0] = point[0] - cp[0];
					vv[1] = point[1] - cp[1];
					vv[2] = point[2] - cp[2];

					double normal[3];
					normals->GetTuple(cid, normal);

					td = vtkMath::Dot(vv, normal);
					distances.push_back(td);
					distIDS.push_back(cid);
				}

				// HACK! We are using the ID of the median point in the lookup, but the distance is the average
				int medIndx = NumberOfDistances / 2;
				CGeneralUtils::Sort2Vectors(distances, distIDS);
			
				double meandist = 0;
				double sdevdist = 0;
				CGeneralUtils::MeanAndSdev(distances, meandist, sdevdist);
				int closestID = distIDS[medIndx];

				newScalars->SetValue(offset, -meandist);
				neighPts->Delete();
				if (CreateReferenceVolume)
				{
					ReferenceVolumeData->SetValue(offset, closestID);
				}
			}
		}
	}
	locator->Delete();
	if (CreateReferenceVolume)
	{
		ReferenceVolumeData->Modified();
	}
}


int vtkOrientedPointSetDistanceFilter::RequestData(
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

	//// need to know the bounding rectangle
	//double bounds[6];
	//for(i=0;i<3;i++)
	//{
	//	bounds[i*2]=input->GetBounds()[i*2];
	//	bounds[i*2+1]=input->GetBounds()[i*2+1];
	//}

	//// estimate the spacing if required
	//if(this->SampleSpacing<=0.0)
	//{
	//	// spacing guessed as cube root of (volume divided by number of points)
	//	this->SampleSpacing = pow((double)(bounds[1]-bounds[0])*
	//		(bounds[3]-bounds[2])*(bounds[5]-bounds[4]) /
	//		(double)NPoints, (double)(1.0/3.0));


	//	this->SampleSpacing /= SampleFactor;

	//	vtkDebugMacro(<<"Estimated sample spacing as: " << this->SampleSpacing );

	//	std::cout <<"Estimated sample spacing as: " << this->SampleSpacing << std::endl;
	//}

	//// allow a border around the volume to allow sampling around the extremes
	//for(i=0;i<3;i++)
	//{
	//	bounds[i*2]-=this->SampleSpacing*3;
	//	bounds[i*2+1]+=this->SampleSpacing*3;
	//}

	//double topleft[3] = {bounds[0],bounds[2],bounds[4]};
	//double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
	//int dim[3];
	//for(i=0;i<3;i++)
	//{
	//	dim[i] = (int)((bottomright[i]-topleft[i])/this->SampleSpacing);
	//}

	//std::cout <<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ") size: " << dim[0] * dim[1] * dim[2]  << std::endl;

	//vtkDebugMacro(<<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" );

	//// initialise the output volume
	//this->GetOutput()->SetWholeExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	//this->GetOutput()->SetUpdateExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);

	//vtkImageData *output = this->AllocateOutputData(outp);

	//vtkDoubleArray *newScalars = 
	//	vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	//output->SetSpacing(this->SampleSpacing, this->SampleSpacing,
	//	this->SampleSpacing);

	//output->SetOrigin(topleft);

	const int NPoints = input->GetNumberOfPoints();

	// need to know the bounding rectangle
	double bounds[6];
	for(int i = 0; i < 3; i++)
	{
		bounds[i*2]=input->GetBounds()[i*2];
		bounds[i*2+1]=input->GetBounds()[i*2+1];
	}

	// estimate the spacing if required
	if(this->SampleSpacing<=0.0)
	{
		// spacing guessed as cube root of (volume divided by number of points)
		this->SampleSpacing = pow(static_cast<double>(bounds[1]-bounds[0])*
			(bounds[3]-bounds[2])*(bounds[5]-bounds[4]) /
			static_cast<double>(NPoints),
			static_cast<double>(1.0/3.0));
		std::cout <<"Estimated sample spacing as: " << this->SampleSpacing << std::endl;
		vtkDebugMacro(<<"Estimated sample spacing as: " << this->SampleSpacing );
	}

	// allow a border around the volume to allow sampling around the extremes
	for(int i = 0; i < 3; i++)
	{
		bounds[i*2]-=this->SampleSpacing * 3;
		bounds[i*2+1]+=this->SampleSpacing * 3;
	}

	double topleft[3] = {bounds[0],bounds[2],bounds[4]};
	double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>((bottomright[i]-topleft[i])/this->SampleSpacing);
	}

	vtkDebugMacro(<<"Created output volume of dimensions: ("
		<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" );

	// initialise the output volume
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	output->SetExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	output->AllocateScalars(outInfo);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),	0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);

	vtkDoubleArray *newScalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
	outInfo->Set(vtkDataObject::SPACING(),	this->SampleSpacing, this->SampleSpacing, this->SampleSpacing);
	outInfo->Set(vtkDataObject::ORIGIN(), topleft, 3);

	vtkDataArray *normals = input->GetPointData()->GetNormals(); 

	if (!normals)
	{
		std::cerr << "No normals found" << std::endl;
	}

	if (CreateReferenceVolume)
	{
		int size = dim[0] * dim[1] * dim[2];
		ReferenceVolumeData = vtkIntArray::New();
		ReferenceVolumeData->SetNumberOfComponents(1);
		ReferenceVolumeData->SetNumberOfValues(size);

		ReferenceVolume = vtkImageData::New();
		ReferenceVolume->CopyStructure(output);
		ReferenceVolume->GetPointData()->SetScalars(ReferenceVolumeData);
	}

	if (DistanceMode == VTK_ORIENTEDDISTANCE_SIMPLE)
		UseSimpleCloseNeighbour(input, dim, newScalars, topleft, normals);	
	else if (DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED)
		UseProjectedCloseNeighbour(input, dim, newScalars, topleft, normals);	
	else if (DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED_MEDIAN)
		UseProjectedCloseNeighbourMedian(input, dim, newScalars, topleft, normals);	
	else if (DistanceMode == VTK_ORIENTEDDISTANCE_PROJECTED_AVERAGE)
		UseProjectedCloseNeighbourAverage(input, dim, newScalars, topleft, normals);	

	return 1;
}

void vtkOrientedPointSetDistanceFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);

	os << indent << "Sample Spacing:" << this->SampleSpacing << "\n";
}

