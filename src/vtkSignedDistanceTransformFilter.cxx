#include "vtkSignedDistanceTransformFilter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkCellLocator.h"
#include "vtkPoints.h"
#include "vtkPlane.h"
#include "vtkTriangle.h"
#include <algorithm>
#include <iostream>
#include <deque>
#include <vtkPolyData.h>

vtkStandardNewMacro(vtkSignedDistanceTransformFilter);

vtkSignedDistanceTransformFilter::vtkSignedDistanceTransformFilter()
{
	this->SampleSpacing = -1.0; // negative values cause the algorithm to make a reasonable guess
	this->SDMMode = VTK_SDM_MODE_LOCALCELL;
	m_MinBoundLength = 1.5;
}


int vtkSignedDistanceTransformFilter::FillInputPortInformation(
  int vtkNotUsed( port ), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}


int vtkSignedDistanceTransformFilter::RequestInformation (
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


vtkDoubleArray * vtkSignedDistanceTransformFilter::AllocateVolume(vtkDataSet *input, vtkImageData *output, vtkInformation *outInfo, double Boundfactor)
{
	// need to know the bounding rectangle
//	double bounds[6];
	for(int i = 0; i < 3; i++)
	{
		bounds[i*2]=input->GetBounds()[i*2];
		bounds[i*2+1]=input->GetBounds()[i*2+1];
	}

	const int NPoints = input->GetNumberOfPoints();

	// estimate the spacing if required
	if(this->SampleSpacing<=0.0)
	{
		// spacing guessed as cube root of (volume divided by number of points)
		this->SampleSpacing = pow(static_cast<double>(bounds[1]-bounds[0])*
			(bounds[3]-bounds[2])*(bounds[5]-bounds[4]) /
			static_cast<double>(NPoints),
			static_cast<double>(1.0/3.0));

		vtkDebugMacro(<<"Estimated sample spacing as: " << this->SampleSpacing );
	}

	// allow a border around the volume to allow sampling around the extremes
	for(int i = 0; i < 3; i++)
	{
		bounds[i*2]-=this->SampleSpacing * Boundfactor;
		bounds[i*2+1]+=this->SampleSpacing * Boundfactor;
		topleft[i] = bounds[i * 2];
		bottomright[i] = bounds[i * 2 + 1];
	}

//	double topleft[3] = {bounds[0],bounds[2],bounds[4]};
//	double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
//	int dim[3];
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
	
	return newScalars;
}



int vtkSignedDistanceTransformFilter::RequestData(
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

  	if (SDMMode == VTK_SDM_MODE_BRUTEFORCE)
		ComputeSDMBruteforce(inputVector, outputVector);
	else if (SDMMode == VTK_SDM_MODE_LOCALCELL)
		ComputeSDMByCellProbing(inputVector, outputVector);
	else if (SDMMode == VTK_SDM_MODE_REALBRUTE)
		ComputeSDMReallyBruteforce(inputVector, outputVector);
	return 1;
}

////-----------------------------------------------------------------------------
//void vtkSignedDistanceTransformFilter::ExecuteData(vtkDataObject *outp)
//{
//	if (SDMMode == VTK_SDM_MODE_BRUTEFORCE)
//		ComputeSDMBruteforce(outp);
//	else if (SDMMode == VTK_SDM_MODE_LOCALCELL)
//		ComputeSDMByCellProbing(outp);
//	else if (SDMMode == VTK_SDM_MODE_REALBRUTE)
//		ComputeSDMReallyBruteforce(outp);
//}


//-----------------------------------------------------------------------------
void vtkSignedDistanceTransformFilter::ComputeSDMByCellProbing(vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	// Initialise the variables we need within this function
//	vtkDataSet *input = this->GetInput();
	
	unsigned int j;
	int i;
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];
	
	  // get the input
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	// get the output
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkDoubleArray *newScalars = AllocateVolume(input, output, outInfo, 10);

	//const int NPoints = input->GetNumberOfPoints();
	//// need to know the bounding rectangle
	//double bounds[6];
	//for(i=0;i<3;i++)
	//{
	//	bounds[i*2]=input->GetBounds()[i*2];
	//	bounds[i*2+1]=input->GetBounds()[i*2+1];
	//}
	//
	//// estimate the spacing if required
	//if(this->SampleSpacing<=0.0)
	//{
	//	// spacing guessed as cube root of (volume divided by number of points)
	//	this->SampleSpacing = pow((double)(bounds[1]-bounds[0])*
	//		(bounds[3]-bounds[2])*(bounds[5]-bounds[4]) /
	//		(double)NPoints, (double)(1.0/3.0));
	//	
	//	vtkDebugMacro(<<"Estimated sample spacing as: " << this->SampleSpacing );
	//}
	//
	//// allow a border around the volume to allow sampling around the extremes
	//for(i=0;i<3;i++)
	//{
	//	bounds[i*2]-=this->SampleSpacing*10;
	//	bounds[i*2+1]+=this->SampleSpacing*10;
	//}
	//
	//double topleft[3] = {bounds[0],bounds[2],bounds[4]};
	//double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
	//int dim[3];
	//for(i=0;i<3;i++)
	//{
	//	dim[i] = (int)((bottomright[i]-topleft[i])/this->SampleSpacing);
	//}
	//
	//std::cout <<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" << std::endl;

	//vtkDebugMacro(<<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" );
	//
	//// initialise the output volume
	//this->GetOutput()->SetWholeExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	//this->GetOutput()->SetUpdateExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	//
	//vtkImageData *output = this->AllocateOutputData(outp);
	//
	//vtkDoubleArray *newScalars = 
	//	vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
	//
	//output->SetSpacing(this->SampleSpacing, this->SampleSpacing,
	//	this->SampleSpacing);
	//
	//output->SetOrigin(topleft);

	double dumMax = sqrt(VTK_DOUBLE_MAX);
	
	std::cout << "Clearing data" << std::endl;
	// Clear data and set them to a dummy value
	for (i = 0; i < newScalars->GetNumberOfTuples(); i++)
	{
		newScalars->SetValue(i, dumMax);
	}

	std::cout << "Marking narrow band with " << input->GetNumberOfCells() << " cells." << std::endl;

	vtkDataArray *normals = input->GetPointData()->GetNormals(); 

	// go through all cells in shape and mark narrow band
	for (i = 0; i < input->GetNumberOfCells(); i++)
	{
		if ((i % 500) == 0)
		{
			std::cout << i << " : ";
		}

		vtkCell *cell = input->GetCell(i);
/*
		{
			int id0 = cell->GetPointId(0);
			int id1 = cell->GetPointId(1);
			int id2 = cell->GetPointId(2);
			if (id0 == 3976 || id1 == 3976 ||  id2 == 3976)
			{
				std::cout << "blingbling" << std::endl;
			}
		}
*/

		double bounds[6];
		cell->GetBounds(bounds);
		
		// Resize bounds to double size
		double minx = bounds[0]; 
		double maxx = bounds[1];
		double lx = maxx-minx;
		double midx = (maxx+minx) / 2;
		double miny = bounds[2]; 
		double maxy = bounds[3];
		double ly = maxy-miny;
		double midy = (maxy+miny) / 2;
		double minz = bounds[4]; 
		double maxz = bounds[5];
		double lz = maxz-minz;
		double midz = (maxz+minz) / 2;

		const double upfactor = 1.5;
		lx = std::max(m_MinBoundLength, lx * upfactor);
		ly = std::max(m_MinBoundLength, ly * upfactor);
		lz = std::max(m_MinBoundLength, lz * upfactor);

//		double lmax = std::max(lx, std::max(lz,ly)) * 2;

		bounds[0] = midx - lx;
		bounds[1] = midx + lx;
		bounds[2] = midy - ly;
		bounds[3] = midy + ly;
		bounds[4] = midz - lz;
		bounds[5] = midz + lz;

		// Find cells in bounds
		double minxt = (bounds[0] - topleft[0]) / SampleSpacing;
		int xl = std::max((int)minxt, 0);
		double maxxt = (bounds[1] - topleft[0]) / SampleSpacing;
		int xh = std::min((int)maxxt+1, dim[0]-1);

		double minyt = (bounds[2] - topleft[1]) / SampleSpacing;
		int yl = std::max((int)minyt, 0);
		double maxyt = (bounds[3] - topleft[1]) / SampleSpacing;
		int yh = std::min((int)maxyt+1, dim[1]-1);

		double minzt = (bounds[4] - topleft[2]) / SampleSpacing;
		int zl = std::max((int)minzt, 0);
		double maxzt = (bounds[5] - topleft[2]) / SampleSpacing;
		int zh = std::min((int)maxzt+1, dim[2]-1);


		for(z = zl; z < zh; z++)
		{
			zOffset = z*dim[1]*dim[0];

			for(y = yl; y < yh; y++)
			{
				yOffset = y*dim[0] + zOffset;
				
				for(x = xl; x < xh; x++)
				{
					offset = x + yOffset;

					point[0] = topleft[0] + x*this->SampleSpacing;
					point[1] = topleft[1] + y*this->SampleSpacing;
					point[2] = topleft[2] + z*this->SampleSpacing;

					// check if points are outside bounding planes
					bool outside = false;
					// First check against bounding planes
					for (j = 0; j < m_BoundingPlanes.size() && !outside; j++)
					{
						vtkPlane *plane = m_BoundingPlanes[j];

						if (plane->EvaluateFunction(point) <= 0)
						{
							outside = true;
						}
					}
					
					if (!outside)
					{
						vtkTriangle *tri = vtkTriangle::SafeDownCast(cell);
						
						if (tri)
						{
							// Find distance from point to current cell
							double pcoords[3];
							double cp[3];
							double dist2;
							double weights[3];
							int subid;
							int res = tri->EvaluatePosition(point, cp, subid, pcoords, dist2, weights);
//							if (res == 1)
//							{
//								std::cout << "res == 1" << std::endl;
//							}
//							// \todo Testing of comparing with m_MinBoundLength
//							if (res != -1 && dist2 < m_MinBoundLength * 2)
							if (res != -1)
							{
								double olddist2 = newScalars->GetValue(offset);
								olddist2 *= olddist2;

								// We only compare the numerical value
								// only update if dist is smaller
								if (dist2 < olddist2)
								{
									// Find vector from point on cell to actual point
									double vv[3];
									vv[0] = point[0] - cp[0];
									vv[1] = point[1] - cp[1];
									vv[2] = point[2] - cp[2];

									double factor = 1;

									// Take closest normal
									int largest = 0;
									if (weights[1] > weights[0])
										largest = 1;
									if (weights[2] > weights[largest])
										largest = 2;

									int id0 = tri->GetPointId(largest);

									double normal[3];
									normals->GetTuple(id0, normal);

/* overkill to use normal interpolation.
									int id0 = tri->GetPointId(0);
									int id1 = tri->GetPointId(1);
									int id2 = tri->GetPointId(2);

									double norm0[3];
									normals->GetTuple(id0, norm0);
									double norm1[3];
									normals->GetTuple(id1, norm1);
									double norm2[3];
									normals->GetTuple(id2, norm2);

									double normal[3];
									// Find x coordinate by weighting x coordinates off the three normals (doing the same with y and z)
									normal[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
									normal[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
									normal[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];
*/
									factor = vtkMath::Dot(vv, normal);

									factor /= fabs(factor);
									newScalars->SetValue(offset, factor * sqrt(dist2));

//									if (x == 73 && y == 60 && z == 97)
//									{
//										std::cout << "dist: " << sqrt(dist2) << " factor: " << factor << 
//											" Cell " << i << " cp: " << cp[0] << " " << cp[1] << " " << cp[2] << " " 
//											<< " P: " << point[0] << " " << point[1] << " " << point[2] << " " << std::endl;
//									}
								}
							}
						}
						else
						{
							std::cerr << "Only triangles are supported" << std::endl;
						}
					}
				}
			}
		}
	}
	std::cout << std::endl;
}

//-----------------------------------------------------------------------------
void vtkSignedDistanceTransformFilter::ComputeSDMBruteforce(vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	// Initialise the variables we need within this function
//	vtkDataSet *input = this->GetInput();
	
	int i;
	unsigned int j;
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];
	
	  // get the input
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	// get the output
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


	//const int NPoints = input->GetNumberOfPoints();
	//// need to know the bounding rectangle
	//double bounds[6];
	//for(i=0;i<3;i++)
	//{
	//	bounds[i*2]=input->GetBounds()[i*2];
	//	bounds[i*2+1]=input->GetBounds()[i*2+1];
	//}
	//
	//// estimate the spacing if required
	//if(this->SampleSpacing<=0.0)
	//{
	//	// spacing guessed as cube root of (volume divided by number of points)
	//	this->SampleSpacing = pow((double)(bounds[1]-bounds[0])*
	//		(bounds[3]-bounds[2])*(bounds[5]-bounds[4]) /
	//		(double)NPoints, (double)(1.0/3.0));
	//	
	//	vtkDebugMacro(<<"Estimated sample spacing as: " << this->SampleSpacing );
	//}
	//
	//// allow a border around the volume to allow sampling around the extremes
	//for(i=0;i<3;i++)
	//{
	//	bounds[i*2]-=this->SampleSpacing*10;
	//	bounds[i*2+1]+=this->SampleSpacing*10;
	//}
	//
	//double topleft[3] = {bounds[0],bounds[2],bounds[4]};
	//double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
	//int dim[3];
	//for(i=0;i<3;i++)
	//{
	//	dim[i] = (int)((bottomright[i]-topleft[i])/this->SampleSpacing);
	//}
	//
	//vtkDebugMacro(<<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" );
	//
	//// initialise the output volume
	//this->GetOutput()->SetWholeExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	//this->GetOutput()->SetUpdateExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	//
	//vtkImageData *output = this->AllocateOutputData(outp);
	//
	//vtkDoubleArray *newScalars = 
	//	vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
	//
	//output->SetSpacing(this->SampleSpacing, this->SampleSpacing,
	//	this->SampleSpacing);
	//
	//output->SetOrigin(topleft);
	
	vtkDoubleArray *newScalars = AllocateVolume(input, output, outInfo, 10);

	vtkCellLocator *locator = vtkCellLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();


	std::cout << "Clearing data" << std::endl;
	// Clear data and set them to a dummy value
	for (i = 0; i < newScalars->GetNumberOfTuples(); i++)
	{
		newScalars->SetValue(i, 100000);
	}

	std::cout << "Marking narrow band with " << input->GetNumberOfCells() << " cells." << std::endl;

	// go through all cells in shape and mark narrow band
	for (i = 0; i < input->GetNumberOfCells(); i++)
	{
		vtkCell *cell = input->GetCell(i);
		double bounds[6];
		cell->GetBounds(bounds);
		
		// Resize bounds to double size
		double minx = bounds[0]; 
		double maxx = bounds[1];
		double lx = maxx-minx;
		double midx = (maxx+minx) / 2;
		double miny = bounds[2]; 
		double maxy = bounds[3];
		double ly = maxy-miny;
		double midy = (maxy+miny) / 2;
		double minz = bounds[4]; 
		double maxz = bounds[5];
		double lz = maxz-minz;
		double midz = (maxz+minz) / 2;
		double lmax = std::max(lx, std::max(lz,ly)) * 2;

		bounds[0] = midx - lmax;
		bounds[1] = midx + lmax;
		bounds[2] = midy - lmax;
		bounds[3] = midy + lmax;
		bounds[4] = midz - lmax;
		bounds[5] = midz + lmax;

		// Find cells in bounds
		double minxt = (bounds[0] - topleft[0]) / SampleSpacing;
		int xl = std::max((int)minxt, 0);
		double maxxt = (bounds[1] - topleft[0]) / SampleSpacing;
		int xh = std::min((int)maxxt+1, dim[0]-1);
		double minyt = (bounds[2] - topleft[1]) / SampleSpacing;
		int yl = std::max((int)minyt, 0);
		double maxyt = (bounds[3] - topleft[1]) / SampleSpacing;
		int yh = std::min((int)maxyt+1, dim[1]-1);
		double minzt = (bounds[4] - topleft[2]) / SampleSpacing;
		int zl = std::max((int)minzt, 0);
		double maxzt = (bounds[5] - topleft[2]) / SampleSpacing;
		int zh = std::min((int)maxzt+1, dim[2]-1);

		for(z = zl; z < zh; z++)
		{
			zOffset = z*dim[1]*dim[0];

			for(y = yl; y < yh; y++)
			{
				yOffset = y*dim[0] + zOffset;
				
				for(x = xl; x < xh; x++)
				{
					offset = x + yOffset;

					point[0] = topleft[0] + x*this->SampleSpacing;
					point[1] = topleft[1] + y*this->SampleSpacing;
					point[2] = topleft[2] + z*this->SampleSpacing;

					// check if points are outside bounding planes
					bool outside = false;
					// First check against bounding planes
					for (j = 0; j < m_BoundingPlanes.size() && !outside; j++)
					{
						vtkPlane *plane = m_BoundingPlanes[j];

						if (plane->EvaluateFunction(point) <= 0)
						{
							outside = true;
						}
					}
					
					if (!outside)
					{
						// Set magic number (== 1)
						newScalars->SetValue(offset, 1);
					}
				}
			}
		}
	}

	// Keep track of seed points for flood fill
	std::deque<XYZPoint> pointque;

	int NarrowNumber = 0;
	std::cout << "Probing values in narrow band" << std::endl;
	// go through the array probing the values
	for(z=0;z<dim[2];z++)
	{
		std::cout << "z: " << z << "/" << dim[2] << " ";
		
		zOffset = z*dim[1]*dim[0];
		point[2] = topleft[2] + z*this->SampleSpacing;
		
		for(y=0;y<dim[1];y++)
		{
			yOffset = y*dim[0] + zOffset;
			point[1] = topleft[1] + y*this->SampleSpacing;
			
			for(x=0;x<dim[0];x++)
			{
				offset = x + yOffset;

				// Check if point is in narrow band and distance shall be calculated
				if (newScalars->GetValue(offset) == 1)
				{
					NarrowNumber++;

					pointque.push_back(XYZPoint(x, y, z));

					point[0] = topleft[0] + x*this->SampleSpacing;
					
					// closest point
					double cp[3];
					
					int sub_id;
					double dist2 = 0;
					vtkIdType cell_id;
					
					locator->FindClosestPoint(point, cp, cell_id, sub_id, dist2);
					
					double td = sqrt(dist2);
					
					// Find vector from point on surface to actual point
					double vv[3];
					vv[0] = point[0] - cp[0];
					vv[1] = point[1] - cp[1];
					vv[2] = point[2] - cp[2];

					double factor = 1;

					// Find weighted normal
					vtkTriangle *tri = vtkTriangle::SafeDownCast(input->GetCell(cell_id));
					if (tri)
					{
						// Closest point in evaluation
						double p2[3];

						// Parametric coordinates
						double pcoords[3];

						// Distance
						double dist2;

						double weights[3];

						int inside = tri->EvaluatePosition(point,p2,sub_id,pcoords,dist2,weights);

						int id0 = tri->GetPointId(0);
						int id1 = tri->GetPointId(1);
						int id2 = tri->GetPointId(2);

						vtkDataArray *normals = input->GetPointData()->GetNormals(); 
						double norm0[3];
						normals->GetTuple(id0, norm0);
						double norm1[3];
						normals->GetTuple(id1, norm1);
						double norm2[3];
						normals->GetTuple(id2, norm2);

						double normal[3];
						// Find x coordinate by weighting x coordinates off the three normals (doing the same with y and z)
						normal[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
						normal[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
						normal[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];

						factor = vtkMath::Dot(vv, normal);

						// Normalise
						factor /= fabs(factor);
					}
					newScalars->SetValue(offset, factor*td);
				}
			}
		}
	}

	std::cout << "Probed " << NarrowNumber << " values in narrow band" << std::endl;
	std::cout << "Compared to " << dim[0] * dim[1] * dim[2] << " values in total volume" << std::endl;

	locator->Delete();

	std::cout << "Flood filling the rest" << std::endl;

	do {
		
		XYZPoint pt = pointque.front();
		pointque.pop_front();
		offset =  pt.x +  pt.y * dim[0] + pt.z * dim[1]*dim[0];
		double vt = newScalars->GetValue(offset);

		// Check neighbours
		for (int xt = -1; xt <= 1; xt++)
		{
			for (int yt = -1; yt <= 1; yt++)
			{
				for (int zt = -1; zt <= 1; zt++)
				{
					XYZPoint cp(pt.x + xt, pt.y + yt, pt.z + zt);
					cp.x = std::max(0, std::min(cp.x, dim[0]-1));
					cp.y = std::max(0, std::min(cp.y, dim[1]-1));
					cp.z = std::max(0, std::min(cp.z, dim[2]-1));

					point[0] = topleft[0] + cp.x*this->SampleSpacing;
					point[1] = topleft[1] + cp.y*this->SampleSpacing;
					point[2] = topleft[2] + cp.z*this->SampleSpacing;

					// check if points are outside bounding planes
					bool outside = false;
					// First check against bounding planes
					for (unsigned int j = 0; j < m_BoundingPlanes.size() && !outside; j++)
					{
						vtkPlane *plane = m_BoundingPlanes[j];

						if (plane->EvaluateFunction(point) <= 0)
						{
							outside = true;
						}
					}

					if (!outside)
					{
						offset =  cp.x +  cp.y * dim[0] + cp.z * dim[1]*dim[0];

						if (newScalars->GetValue(offset) > 99999)
						{
							newScalars->SetValue(offset, vt);
							pointque.push_back(cp);
						}
					}
				}
			}
		}
	} while(pointque.size());
	std::cout << "Flood fill finished" << std::endl;

}


//-----------------------------------------------------------------------------
void vtkSignedDistanceTransformFilter::ComputeSDMReallyBruteforce(vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	std::cout << "Filling distance map brute force" << std::endl;

	// Initialise the variables we need within this function
//	vtkDataSet *input = this->GetInput();

//	int i;
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	  // get the input
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	// get the output
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkDataArray *normals = input->GetPointData()->GetNormals(); 
	if (!normals)
	{
		vtkErrorMacro(<<"Normals not defined");
		return;
	}

	//const int NPoints = input->GetNumberOfPoints();
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

	//	vtkDebugMacro(<<"Estimated sample spacing as: " << this->SampleSpacing );

	//	std::cout <<"Estimated sample spacing as: " << this->SampleSpacing << std::endl;
	//}

	//// allow a border around the volume to allow sampling around the extremes
	//for(i=0;i<3;i++)
	//{
	//	//bounds[i*2]-=this->SampleSpacing*10;
	//	//bounds[i*2+1]+=this->SampleSpacing*10;
	//	bounds[i*2]-=this->SampleSpacing*2;
	//	bounds[i*2+1]+=this->SampleSpacing*2;
	//}

	//double topleft[3] = {bounds[0],bounds[2],bounds[4]};
	//double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
	//int dim[3];
	//for(i=0;i<3;i++)
	//{
	//	dim[i] = (int)((bottomright[i]-topleft[i])/this->SampleSpacing);
	//}

	//vtkDebugMacro(<<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" );

	//std::cout <<"Created output volume of dimensions: ("
	//	<< dim[0] << ", " << dim[1] << ", " << dim[2] << ")" << std::endl;


	//// initialise the output volume
	//this->GetOutput()->SetWholeExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
	//this->GetOutput()->SetUpdateExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);

	//vtkImageData *output = this->AllocateOutputData(outp);

	//vtkDoubleArray *newScalars = 
	//	vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	//output->SetSpacing(this->SampleSpacing, this->SampleSpacing,
	//	this->SampleSpacing);

	//output->SetOrigin(topleft);

	vtkDoubleArray *newScalars = AllocateVolume(input, output, outInfo, 2);

	vtkCellLocator *locator = vtkCellLocator::New();
	locator->SetDataSet(input);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();

	vtkIdList *neighbors = vtkIdList::New();
	neighbors->Allocate(VTK_CELL_SIZE);

	vtkPolyData *pd = vtkPolyData::SafeDownCast(input);
	if (!pd)
	{
		std::cerr << "Could not transform input to polydata" << std::endl;
		return;
	}
	pd->BuildLinks();

	// go through the array probing the values
	for(z=0;z<dim[2];z++)
	{
		std::cout << "z: " << z << "/" << dim[2] << " ";

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

				int sub_id;
				double dist2 = 0;
				vtkIdType cell_id;

				locator->FindClosestPoint(point, cp, cell_id, sub_id, dist2);

				double td = sqrt(dist2);

				// Find vector from point on surface to actual point
				double vv[3];
				vv[0] = point[0] - cp[0];
				vv[1] = point[1] - cp[1];
				vv[2] = point[2] - cp[2];

				double factor = 1;

				// Find weighted normal
				vtkTriangle *tri = vtkTriangle::SafeDownCast(input->GetCell(cell_id));
				if (tri)
				{
					int id0 = tri->GetPointId(0);
					int id1 = tri->GetPointId(1);
					int id2 = tri->GetPointId(2);

					// First check if it an edge triangle
					pd->GetCellEdgeNeighbors(cell_id,id0,id1, neighbors);
					int numNei = neighbors->GetNumberOfIds();
					if (numNei == 1)
					{
						pd->GetCellEdgeNeighbors(cell_id,id1,id2, neighbors);
						numNei = neighbors->GetNumberOfIds();
						if (numNei == 1)
						{
							pd->GetCellEdgeNeighbors(cell_id,id0,id2, neighbors);
							numNei = neighbors->GetNumberOfIds();
							if (numNei == 1)
							{
								// Closest point in evaluation
								double p2[3];

								// Parametric coordinates
								double pcoords[3];

								// Distance
								double dist2;

								double weights[3];

								int inside = tri->EvaluatePosition(point,p2,sub_id,pcoords,dist2,weights);

								vtkDataArray *normals = input->GetPointData()->GetNormals(); 
								double norm0[3];
								normals->GetTuple(id0, norm0);
								double norm1[3];
								normals->GetTuple(id1, norm1);
								double norm2[3];
								normals->GetTuple(id2, norm2);

								double normal[3];
								// Find x coordinate by weighting x coordinates off the three normals (doing the same with y and z)
								normal[0] = weights[0] * norm0[0] + weights[1] * norm1[0] + weights[2] * norm2[0];
								normal[1] = weights[0] * norm0[1] + weights[1] * norm1[1] + weights[2] * norm2[1];
								normal[2] = weights[0] * norm0[2] + weights[1] * norm1[2] + weights[2] * norm2[2];

								factor = vtkMath::Dot(vv, normal);

								// Normalise
								factor /= fabs(factor);
							}
						}
					}
				}
				newScalars->SetValue(offset, -factor*td);
			}
		}
	}
	locator->Delete();
	neighbors->Delete();
}

void vtkSignedDistanceTransformFilter::InsertBoundingplane(vtkPlane *p)
{
	m_BoundingPlanes.push_back(p);
}

void vtkSignedDistanceTransformFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
	
	os << indent << "Sample Spacing:" << this->SampleSpacing << "\n";
}
