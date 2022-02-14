#include "vtkDistanceFieldMRFRegularisation.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkpointdata.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <algorithm>
#include <vnl/algo/vnl_conjugate_gradient.h>
#include <vnl/algo/vnl_lbfgs.h>
#include <thread>
#include <future>

#ifdef USECHOLMOD
#include "cholmod.h"
#endif
#include "GeneralUtils.h"

vtkStandardNewMacro(vtkDistanceFieldMRFRegularisation);

vtkDistanceFieldMRFRegularisation::vtkDistanceFieldMRFRegularisation()
{
	WeightValueHigh = 3;
	Iterations = 1;
	OrgScalars = vtkDoubleArray::New();
	OrgScalars->SetNumberOfComponents(1);
	GlobalBeta = 0.5;
	ReferenceVolume = NULL;
	ReferencePD = NULL;
	LocalWeights = NULL;
	UseLocalWeighting = false;
	MinRMS = 1e-6;
	WeightVol = NULL;
	WeightType = VTK_WEIGHT_BASED_ON_INPUT;
//	ChangeVolume = NULL;
	PredictionStep = 10;
	PriorType = VTK_PRIOR_ENERGY_DIFF_LAPLACE;
	DeleteWeightVol = true;
	CGTolerance = 10;
	Optimisation = VTK_OPTIMISATION_ICM;
	BandedICM = 0;
	PerVoxelICMMemory = 0;
	ICMUpdateVolume = NULL;
	NumThreads = 0;
	m_ConvolutionCoefficients = NULL;
	FinalSplineEnergy = 0;
}

vtkDistanceFieldMRFRegularisation::~vtkDistanceFieldMRFRegularisation()
{
	if (OrgScalars)
		OrgScalars->Delete();
	if (DeleteWeightVol && LocalWeights)
		LocalWeights->Delete();
	if (DeleteWeightVol && WeightVol)
		WeightVol->Delete();
// 	if (ChangeVolume)
// 		ChangeVolume->Delete();
	if (ICMUpdateVolume)
		ICMUpdateVolume->Delete();
	if (m_ConvolutionCoefficients)
		delete [] m_ConvolutionCoefficients;
	
}

void vtkDistanceFieldMRFRegularisation::VisualiseCoefficients()
{
	std::string outputname = "C:/rrplocal/data/IMM/Ohtake/MPUTest/DoubleLaplaceFilter.vtk";

	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);

	for (unsigned int i = 0; i < m_FilterCoefficients.size(); i++)
	{
		CFilterCoefficient c = m_FilterCoefficients[i];

		double p[3];
		p[0] = c.p[0]; p[1] = c.p[1]; p[2] = c.p[2]; 

		int id = pts->InsertNextPoint(p);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);
		scalars->InsertNextTuple(&c.v);
	}
	vtkPolyData *pd = vtkPolyData::New();
	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->SetScalars(scalars);
	pts->Delete();
	verts->Delete();
	scalars->Delete();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInputData(pd);
	writer->SetFileName(outputname.c_str());
	writer->Write();

	writer->Delete();
	pd->Delete();
}

void vtkDistanceFieldMRFRegularisation::CreateCoefficients()
{
	m_FilterCoefficients.clear();
	
	double scale = -m_CoeffCube[2][2][2];
	if (scale == 0)
	{
		std::cerr << "Something is wrong" << std::endl;
		return;
	}

	// Create convolution kernel
	m_ConvolutionCoefficients = new double[125];
	for (int x = 0; x < 5; x++)
	{
		for (int y = 0; y < 5; y++)
		{
			for (int z = 0; z < 5; z++)
			{
				if (m_CoeffCube[x][y][z] != 0 && (x!=2 || y!= 2 || z!=2))
				{
					CFilterCoefficient c;
					c.p[0] = x-2;
					c.p[1] = y-2;
					c.p[2] = z-2;
					c.v = m_CoeffCube[x][y][z] / scale;
					m_FilterCoefficients.push_back(c);
					m_ConvolutionCoefficients[x + y*5 + z * 5 * 5] = m_CoeffCube[x][y][z] / scale;
				}
				else
				{
					m_ConvolutionCoefficients[x + y*5 + z * 5 * 5] = 0;
				}
			}
		}
	}
	std::cout << "Found " << m_FilterCoefficients.size() << " coefficients" << std::endl;


	//for (unsigned int i = 0; i < m_FilterCoefficients.size(); i++)
	//{
	//	CFilterCoefficient c = m_FilterCoefficients[i];
	//	std::cout << c.p[0] << ", " << c.p[1] << ", " << c.p[2] << ", " << std::endl;
	//}
//	VisualiseCoefficients();
}


void vtkDistanceFieldMRFRegularisation::LocalLaplaceCoefficients(int x, int y, int z, double fac)
{
	// Move to center of cube
	int xt = x + 2;
	int yt = y + 2;
	int zt = z + 2;

	// -6 in the center. 1 in the 6 neighbours
	m_CoeffCube[xt][yt][zt] -= 6*fac;
	m_CoeffCube[xt+1][yt][zt] += fac;
	m_CoeffCube[xt-1][yt][zt] += fac;
	m_CoeffCube[xt][yt+1][zt] += fac;
	m_CoeffCube[xt][yt-1][zt] += fac;
	m_CoeffCube[xt][yt][zt+1] += fac;
	m_CoeffCube[xt][yt][zt-1] += fac;
}

void vtkDistanceFieldMRFRegularisation::CreateLaplaceCoefficients()
{
	m_CoeffCube.resize(5);
	for (int x = 0; x < 5; x++)
	{
		m_CoeffCube[x].resize(5);
		for (int y = 0; y < 5; y++)
		{
			m_CoeffCube[x][y].resize(5);
			for (int z = 0; z < 5; z++)
			{
				m_CoeffCube[x][y][z] = 0;
			}
		}
	}

	LocalLaplaceCoefficients( 0, 0, 0, -6.0);
	LocalLaplaceCoefficients(-1, 0, 0, 1.0);
	LocalLaplaceCoefficients( 1, 0, 0, 1.0);
	LocalLaplaceCoefficients( 0, 1, 0, 1.0);
	LocalLaplaceCoefficients( 0,-1, 0, 1.0);
	LocalLaplaceCoefficients( 0, 0, 1, 1.0);
	LocalLaplaceCoefficients( 0, 0,-1, 1.0);

	CreateCoefficients();
}

void vtkDistanceFieldMRFRegularisation::CreateBandedVisitOrder( int dim[3], vtkImageData *infield )
{
	double spacing[3];
	infield->GetSpacing(spacing);
	double SS = spacing[0];

	vtkDoubleArray *InScalars = 
		vtkDoubleArray::SafeDownCast(infield->GetPointData()->GetScalars());

	double maxDist = 6 * SS;
	double minDist = - maxDist;

	m_VisitOrder.clear();

	int x, y, z;
	int zOffset,yOffset,offset;

	for(z = 2; z < dim[2]-2; z++)
	{
		zOffset = z*dim[1]*dim[0];
		for(y = 2;y < dim[1]-2; y++)
		{
			yOffset = y*dim[0] + zOffset;

			for(x = 2;x < dim[0]-2; x++)
			{
				offset = x + yOffset;
				
				double val = InScalars->GetValue(offset);
				if (val < maxDist && val > minDist)
				{
					CVoxelID vid;
					vid.id[0] = x;
					vid.id[1] = y;
					vid.id[2] = z;

					m_VisitOrder.push_back(vid);
				}
			}
		}
	}

	int size = m_VisitOrder.size();

	// Now create random visit order
	for (int i = 0; i < size; i++)
	{
		int swap = vtkMath::Random(0, size);
		CVoxelID t = m_VisitOrder[i];
		m_VisitOrder[i] = m_VisitOrder[swap];
		m_VisitOrder[swap] = t;
	}

	std::cout << "Created swapped visit order vector for ICM with size " << size << " max size was " << dim[0] * dim[1] * dim[2] << std::endl;
}


void vtkDistanceFieldMRFRegularisation::CreateVisitOrder(int dim[3])
{
	int size = dim[0]*dim[1]*dim[2];

	m_VisitOrder.resize(size);

	int x, y, z;
	int zOffset,yOffset,offset;

	for(z = 0; z < dim[2]; z++)
	{
		zOffset = z*dim[1]*dim[0];
		for(y = 0;y < dim[1]; y++)
		{
			yOffset = y*dim[0] + zOffset;

			for(x = 0;x < dim[0]; x++)
			{
				offset = x + yOffset;
				CVoxelID vid;
				vid.id[0] = x;
				vid.id[1] = y;
				vid.id[2] = z;
				
				m_VisitOrder[offset] = vid;
			}
		}
	}


	//unsigned int i;
	//for (i = 0; i < size; i++)
	//{
	//	m_VisitOrder[i] = i;
	//}

	// Now create random visit order
	for (int i = 0; i < size; i++)
	{
		int swap = vtkMath::Random(0, size);
		CVoxelID t = m_VisitOrder[i];
		m_VisitOrder[i] = m_VisitOrder[swap];
		m_VisitOrder[swap] = t;
	}
}


double vtkDistanceFieldMRFRegularisation::GetValue(vtkDoubleArray* data, int dims[3], int x, int y, int z)
{
	// Border conditions
	if (x < 0)
		x = 0;
	if (x > dims[0]-1)
		x = dims[0]-1;
	if (y < 0)
		y = 0;
	if (y > dims[1]-1)
		y = dims[1]-1;
	if (z < 0)
		z = 0;
	if (z > dims[2]-1)
		z = dims[2]-1;

	int offset = z * dims[0] * dims[1] + y * dims[0] + x;
	return data->GetValue(offset);
}

void vtkDistanceFieldMRFRegularisation::SetValue(vtkDoubleArray* data, int dims[3], int x, int y, int z, double val)
{
	int offset = z * dims[0] * dims[1] + y * dims[0] + x;
	data->SetValue(offset, val);
}

double vtkDistanceFieldMRFRegularisation::LocalBandedDoubleLaplacianFiltering( vtkDoubleArray* data, int dims[3], CVoxelID idx )
{
	int x = idx.id[0];
	int y = idx.id[1];
	int z = idx.id[2];

	double val = 0;

	for (unsigned int i = 0; i < m_FilterCoefficients.size(); i++)
	{
		CFilterCoefficient c = m_FilterCoefficients[i];
		int xt = c.p[0];
		int yt = c.p[1];
		int zt = c.p[2];

		int offset = (z+zt) * dims[0] * dims[1] + (y+yt) * dims[0] + x + xt ;
		double r = c.v * data->GetValue(offset);
		val += r;
	}

//	double orgvalue = GetValue(OrgScalars, dims, x, y, z);
	double orgvalue = OrgScalars->GetValue(z * dims[0] * dims[1] + y * dims[0] + x);

	double oldvalue = data->GetValue(z * dims[0] * dims[1] + y * dims[0] + x);
//	double oldvalue = GetValue(data, dims, x, y, z);

	double newval = 0;
	double change = 0;

	//Apply weight
	if (!UseLocalWeighting)
	{
		newval = GlobalBeta * orgvalue + (1-GlobalBeta) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
//		double localAlpha = GetValue(LocalWeights, dims, x, y, z) * GlobalBeta;
		double localAlpha = LocalWeights->GetValue(z * dims[0] * dims[1] + y * dims[0] + x) * GlobalBeta;
		newval = localAlpha * orgvalue + (1-localAlpha) * val;
		change = newval-oldvalue;
	}
	SetValue(data, dims, x, y, z, newval);

	return change;
}

double vtkDistanceFieldMRFRegularisation::LocalDoubleLaplacianFilteringWithSeperateOutput( vtkDoubleArray* data, vtkDoubleArray* outdata, int dims[3], CVoxelID idx )
{
	int x = idx.id[0];
	int y = idx.id[1];
	int z = idx.id[2];

	double val = 0;

	for (unsigned int i = 0; i < m_FilterCoefficients.size(); i++)
	{
		CFilterCoefficient c = m_FilterCoefficients[i];
		int xt = c.p[0];
		int yt = c.p[1];
		int zt = c.p[2];

		double r = c.v * GetValue(data, dims, x + xt, y + yt, z + zt);
		val += r;
	}

	double orgvalue = GetValue(OrgScalars, dims, x, y, z);
	double oldvalue = GetValue(data, dims, x, y, z);

	double newval = 0;
	double change = 0;

	//Apply weight
	if (!UseLocalWeighting)
	{
		newval = GlobalBeta * orgvalue + (1-GlobalBeta) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
		double localAlpha = GetValue(LocalWeights, dims, x, y, z) * GlobalBeta;
		newval = localAlpha * orgvalue + (1-localAlpha) * val;
		change = newval-oldvalue;

		// Test regularisation
		newval  = oldvalue + 0.5 * change;
	}
	SetValue(outdata, dims, x, y, z, newval);

	if (change > 1000 || change < -1000)
	{
		std::cerr << "Something wrong with change " << change << std::endl;
	}

	return change;
}



double vtkDistanceFieldMRFRegularisation::LocalDoubleLaplacianFilteringUsingConvolutionAndBorderCheck( vtkDoubleArray* data, vtkDoubleArray* outdata, int dims[3], int xo, int yo, int zo )
{
	double val = 0;

	for(int z = 0; z < 5; z++)
	{
		int iz = zo + z - 2; 
		if (iz < 0)
			iz = 0;
		if (iz > dims[2]-1)
			iz = dims[2]-1;

		for(int y = 0; y < 5; y++)
		{
			int iy = yo + y - 2; 
			if (iy < 0)
				iy = 0;
			if (iy > dims[1]-1)
				iy = dims[1]-1;

			for(int x = 0; x < 5; x++)
			{
				int filtoffset = x + y * 5 + z * 5 * 5;
				int ix = xo + x - 2; 

				if (ix < 0)
					ix = 0;
				if (ix > dims[0]-1)
					ix = dims[0]-1;
				int imgoffset = iz * dims[0] * dims[1] + iy * dims[0] + ix;
				val += data->GetValue(imgoffset) * m_ConvolutionCoefficients[filtoffset];
			}
		}
	}

	double orgvalue = GetValue(OrgScalars, dims, xo, yo, zo);
	double oldvalue = GetValue(data, dims, xo, yo, zo);

	double newval = 0;
	double change = 0;

	//Apply weight
	if (!UseLocalWeighting)
	{
		newval = GlobalBeta * orgvalue + (1-GlobalBeta) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
		double localAlpha = GetValue(LocalWeights, dims, xo, yo, zo) * GlobalBeta;
		newval = localAlpha * orgvalue + (1-localAlpha) * val;
		change = newval-oldvalue;

		// Test regularisation
		newval  = oldvalue + 0.5 * change;
	}
	int offset = zo * dims[0] * dims[1] + yo * dims[0] + xo;
	outdata->SetValue(offset, newval);

	//DebugVolume[offset]++;

	//if (change > 1000 || change < -1000)
	//{
	//	std::cerr << "Something wrong with change " << change << std::endl;
	//}

	return change;

}

double vtkDistanceFieldMRFRegularisation::LocalDoubleLaplacianFilteringUsingConvolution( vtkDoubleArray* data, vtkDoubleArray* outdata, int dims[3], int xo, int yo, int zo )
{
	double val = 0;

	for(int z = 0; z < 5; z++)
	{
		for(int y = 0; y < 5; y++)
		{
			for(int x = 0; x < 5; x++)
			{
				int filtoffset = x + y * 5 + z * 5 * 5;
				int ix = xo + x - 2; 
				int iy = yo + y - 2; 
				int iz = zo + z - 2; 
				int imgoffset = iz * dims[0] * dims[1] + iy * dims[0] + ix;

				val += data->GetValue(imgoffset) * m_ConvolutionCoefficients[filtoffset];
			}
		}
	}

	double orgvalue = GetValue(OrgScalars, dims, xo, yo, zo);
	double oldvalue = GetValue(data, dims, xo, yo, zo);

	double newval = 0;
	double change = 0;

	//Apply weight
	if (!UseLocalWeighting)
	{
		newval = GlobalBeta * orgvalue + (1-GlobalBeta) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
		double localAlpha = GetValue(LocalWeights, dims, xo, yo, zo) * GlobalBeta;
		newval = localAlpha * orgvalue + (1-localAlpha) * val;
		change = newval-oldvalue;

		// Test regularisation
		newval  = oldvalue + 0.5 * change;
	}
	int offset = zo * dims[0] * dims[1] + yo * dims[0] + xo;
	outdata->SetValue(offset, newval);

	//DebugVolume[offset]++;

	//if (change > 1000 || change < -1000)
	//{
	//	std::cerr << "Something wrong with change " << change << std::endl;
	//}

	return change;
}


double vtkDistanceFieldMRFRegularisation::LocalDoubleLaplacianFiltering( vtkDoubleArray* data, int dims[3], CVoxelID idx )
{
	int x = idx.id[0];
	int y = idx.id[1];
	int z = idx.id[2];

	double val = 0;

	for (unsigned int i = 0; i < m_FilterCoefficients.size(); i++)
	{
		CFilterCoefficient c = m_FilterCoefficients[i];
		int xt = c.p[0];
		int yt = c.p[1];
		int zt = c.p[2];

		double r = c.v * GetValue(data, dims, x + xt, y + yt, z + zt);
		val += r;
	}

	double orgvalue = GetValue(OrgScalars, dims, x, y, z);
	double oldvalue = GetValue(data, dims, x, y, z);

	double newval = 0;
	double change = 0;

	//Apply weight
	if (!UseLocalWeighting)
	{
		newval = GlobalBeta * orgvalue + (1-GlobalBeta) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
		double localAlpha = GetValue(LocalWeights, dims, x, y, z) * GlobalBeta;
		newval = localAlpha * orgvalue + (1-localAlpha) * val;
		change = newval-oldvalue;

		//double localAlpha = GetValue(LocalWeights, dims, x, y, z);
		//newval = localAlpha * GlobalBeta * orgvalue + (1-GlobalBeta) * val;
		//SetValue(data, dims, x, y, z, newval);
		//change = oldvalue-newval;
	}
	//if (boost)
	//{
	//	newval = oldvalue + 100 * change;
	//}
	SetValue(data, dims, x, y, z, newval);

	return change;
}


double vtkDistanceFieldMRFRegularisation::LocalUniformGradientFiltering( vtkDoubleArray* data, int dims[3], CVoxelID idx, double spacing )
{
	double s = spacing;

	int x = idx.id[0];
	int y = idx.id[1];
	int z = idx.id[2];

	double orgvalue = GetValue(OrgScalars, dims, x, y, z);
	double oldvalue = GetValue(data, dims, x, y, z);

	double newval = 0;
	double change = 0;

	double dx = GetValue(data, dims, x-1, y  , z   );
	double dy = GetValue(data, dims, x  , y-1, z   );
	double dz = GetValue(data, dims, x  , y  , z-1 );

	double test = 4 * (dx + dy + dz) * (dx + dy + dz) - 12 * (dx * dx + dy * dy + dz * dz - s * s);
	double val = oldvalue;
	if (test >= 0)
	{
		double val1 = (2 * (dx + dy + dz) - sqrt(test)) / 6.0;
		double val2 = (2 * (dx + dy + dz) + sqrt(test)) / 6.0;

		if (abs(val1) < abs (val2))
//		if (abs(val1-oldvalue) < abs (val2-oldvalue))
			val = val1;
		else
			val = val2;
	}

	if (!UseLocalWeighting)
	{
		// t = 1: only original distance.
		// t = 0: pure smoothing.
		double t  = std::min(GlobalBeta, 1.0);

		newval = t * orgvalue + (1-t) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
		double localAlpha = GetValue(LocalWeights, dims, x, y, z);

		// t = 1: only original distance.
		// t = 0: pure smoothing.
		double t  = std::min(GlobalBeta * localAlpha, 1.0);

		newval = t * orgvalue + (1-t) * val;
		change = newval-oldvalue;
	}

	SetValue(data, dims, x, y, z, newval);
	return change;
}


double vtkDistanceFieldMRFRegularisation::LocalSmoothingFiltering( vtkDoubleArray* data, int dims[3], CVoxelID idx )
{
	//int x = idx.id[0];
	//int y = idx.id[1];
	//int z = idx.id[2];

	//double val = GetValue(data, dims, x-1, y  , z   ) +
	//			 GetValue(data, dims, x+1, y  , z   ) +
	//			 GetValue(data, dims, x  , y-1, z   ) +
	//			 GetValue(data, dims, x  , y+1, z   ) +
	//			 GetValue(data, dims, x  , y  , z-1 ) +
	//			 GetValue(data, dims, x  , y  , z+1 );

	//val /= 6.0;
	//
	//SetValue(data, dims, x, y, z, val);

	double N = 6;

	int x = idx.id[0];
	int y = idx.id[1];
	int z = idx.id[2];

	double orgvalue = GetValue(OrgScalars, dims, x, y, z);
	double oldvalue = GetValue(data, dims, x, y, z);

	double newval = 0;
	double change = 0;

	double val = GetValue(data, dims, x-1, y  , z   ) +
		GetValue(data, dims, x+1, y  , z   ) +
		GetValue(data, dims, x  , y-1, z   ) +
		GetValue(data, dims, x  , y+1, z   ) +
		GetValue(data, dims, x  , y  , z-1 ) +
		GetValue(data, dims, x  , y  , z+1 );

	val /= N;

	if (!UseLocalWeighting)
	{
		// t = 1: only original distance.
		// t = 0: pure smoothing.
		double t  = std::min(GlobalBeta, 1.0);

		newval = t * orgvalue + (1-t) * val;
		change = newval-oldvalue;
	}
	else
	{
		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density
		double localAlpha = GetValue(LocalWeights, dims, x, y, z);

		// t = 1: only original distance.
		// t = 0: pure smoothing.
		double t  = std::min(GlobalBeta * localAlpha, 1.0);

		newval = t * orgvalue + (1-t) * val;
		change = newval-oldvalue;
	}

	SetValue(data, dims, x, y, z, newval);
	return change;
}

void vtkDistanceFieldMRFRegularisation::Median(std::vector<double>& x, double fractile, double& median)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	int med = (int)(n * fractile);
	if (med < 0 || med >= n)
	{
		std::cerr << "fractile incorrect" << std::endl;
		return;
	}

	std::sort(x.begin(), x.end());
	median = x[med];
}

void vtkDistanceFieldMRFRegularisation::DistanceBasedWeights( vtkImageData* input, vtkIntArray * RefVolData)
{
	int size = RefVolData->GetNumberOfTuples();

	double bounds[6];
	input->GetBounds(bounds);
	double spacing[3];
	input->GetSpacing(spacing);
	double SampleSpacing = spacing[0];

	double topleft[3] = {bounds[0],bounds[2],bounds[4]};
	double bottomright[3] = {bounds[1],bounds[3],bounds[5]};
	int dim[3];
	input->GetDimensions(dim);

	// go through the array probing the values
	int x,y,z;
	int zOffset,yOffset,offset;
	double point[3];

	for(z=0;z<dim[2];z++)
	{
		std::cout << "\rz: " << z << "/" << dim[2] << " ";
		zOffset = z*dim[1]*dim[0];
		point[2] = topleft[2] + z*SampleSpacing;

		for(y=0;y<dim[1];y++)
		{
			yOffset = y*dim[0] + zOffset;
			point[1] = topleft[1] + y*SampleSpacing;

			for(x=0;x<dim[0];x++)
			{
				offset = x + yOffset;
				point[0] = topleft[0] + x*SampleSpacing;

				int id = RefVolData->GetValue(offset);

				double p[3];
				ReferencePD->GetPoint(id, p);

				double val = sqrt(vtkMath::Distance2BetweenPoints(point, p));

				LocalWeights->SetValue(offset, val);

//				medHelp[offset] = val;
			}
		}
	}
	//// Calculate 5% and 95% fractilese
	//double frac05;
	//double frac95;

	//Median(medHelp, 0.05, frac05);
	//Median(medHelp, 0.95, frac95);

	//std::cout << "Weight statistics: 5% fractile " << frac05 << " 95% fractile " << frac95 << std::endl;

	// alpha = 1: only original distance. Very high local sampling density
	// alpha = 0: pure smoothing. Very low local sampling density

	double low = 0;
	double high = WeightValueHigh;

	// Normalise weights so they are between 0 and 1
	for (int i = 0; i< size; i++)
	{
		double val = LocalWeights->GetValue(i);

		//			double normval = 1.0-std::max(std::min((val-frac05) / (frac95-frac05), 1.0), 0.0);
		double normval = 1.0-std::max(std::min((val-low) / (high-low), 1.0), 0.0);

		LocalWeights->SetValue(i, normval);
	}
}

void vtkDistanceFieldMRFRegularisation::CreateLocalWeights(vtkImageData* input)
{
	vtkIntArray *RefVolData = vtkIntArray::SafeDownCast(ReferenceVolume->GetPointData()->GetScalars());
	if (!RefVolData)
	{
		std::cerr << "Something wrong with reference volume" << std::endl;
		return;
	}

	int size = RefVolData->GetNumberOfTuples();
	LocalWeights = vtkDoubleArray::New();
	LocalWeights->SetNumberOfComponents(1);
	LocalWeights->SetNumberOfValues(size);

	// Create a helper vector to do the statistics...not very efficient. But only done once
	std::vector<double> medHelp(size);

	if (WeightType == VTK_WEIGHT_BASED_ON_INPUT)
	{
		vtkDoubleArray *values = 
			vtkDoubleArray::SafeDownCast(ReferencePD->GetPointData()->GetScalars());
		if (!values)
		{
			std::cerr << "No scalars in reference PD" << std::endl;
			return;
		}

		// Initially fill the weights with the raw lookup values
		for (int i = 0; i< size; i++)
		{
			int id = RefVolData->GetValue(i);

			double val = values->GetValue(id);

			LocalWeights->SetValue(i, val);

			medHelp[i] = val;
		}
		// Calculate 5% and 95% fractilese
		double frac05;
		double frac95;

		Median(medHelp, 0.05, frac05);
		Median(medHelp, 0.95, frac95);

		std::cout << "Weight statistics: 5% fractile " << frac05 << " 95% fractile " << frac95 << std::endl;

		// alpha = 1: only original distance. Very high local sampling density
		// alpha = 0: pure smoothing. Very low local sampling density

		// Normalise weights so they are between 0 and 1
		for (int i = 0; i< size; i++)
		{
			double val = LocalWeights->GetValue(i);

			double normval = std::max(std::min((val-frac05) / (frac95-frac05), 1.0), 0.0);

			LocalWeights->SetValue(i, normval);
		}
	}
	else
	{
		DistanceBasedWeights(input, RefVolData);
	}

	LocalWeights->Modified();
	UseLocalWeighting = true;
	std::cout << std::endl;
}

void vtkDistanceFieldMRFRegularisation::SimpleExecute(vtkImageData* input,
                                                vtkImageData* output)
{
	vtkDoubleArray *InScalars = 
		vtkDoubleArray::SafeDownCast(input->GetPointData()->GetScalars());

	vtkDoubleArray *OutScalars = 
		vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	if (!InScalars || !OutScalars)
	{
		std::cerr << "Something wrong with input or output" << std::endl;
		return;
	}

	// Determine number of threads
	if (NumThreads != -1)
	{
		// If NumThreads = 0 we automatically determine the number of cores. We leave one core free
		if (NumThreads == 0)
		{	
			NumThreads = std::thread::hardware_concurrency() - 1;
		}
	}

	int dims[3];
	input->GetDimensions(dims);

	int size = dims[0]*dims[1]*dims[2];

	OrgScalars->Resize(size);

	// Copy data
	for(int i = 0; i < size; i++)
	{
		OutScalars->SetValue(i, InScalars->GetValue(i));
		OrgScalars->SetValue(i, InScalars->GetValue(i));
	}

	if (PerVoxelICMMemory)
	{
		// Create a bookkeeping volume
		ICMUpdateVolume = vtkIntArray::New();
		ICMUpdateVolume->Resize(size);
		for (int i = 0; i < size; i++)
		{
			ICMUpdateVolume->SetValue(0,0);
		}
	}

	if (BandedICM)
	{
		CreateBandedVisitOrder(dims, input);
	}
	else
	{
		CreateVisitOrder(dims);
	}

	CreateLaplaceCoefficients();


	if (ReferenceVolume != NULL && ReferencePD != NULL)
	{
		CreateLocalWeights(input);
		WeightVol = vtkImageData::New();
		WeightVol->CopyStructure(input);
		WeightVol->GetPointData()->SetScalars(LocalWeights);
	}


	if (UseLocalWeighting)
	{
		std::cout << "Filtering volume using local weights" << std::endl;
	}
	else
	{
		std::cout << "Filtering volume using a global weight of " << GlobalBeta << std::endl;
	}

	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		Optimisation = VTK_OPTIMISATION_ICM;
	}

	if (Optimisation == VTK_OPTIMISATION_ICM)
	{
		if (PerVoxelICMMemory)
		{
			UseICMOptimisationWithMemory(dims, OutScalars, input);
		}
		else
		{
//			NumThreads = -1;
			// If less than two threads we just use single mode
			if (NumThreads <= 1)
			{
				UseICMOptimisation(dims, OutScalars, input);
//				UseICMOptimisationWithConvolution(dims, OutScalars, input);
			}
			else
			{
//				UseICMOptimisationMultithreaded(dims, OutScalars, input);
//				UseICMOptimisationWithSeperateOutputMultithreaded(dims, OutScalars, input);
				UseICMOptimisationWithConvolutionMultithreaded(dims, OutScalars, input);
			}
		}
	}
	else if (Optimisation == VTK_OPTIMISATION_CONJGRAD)
	{
		UseConjugateGradientOptimisation(dims, OutScalars);
	}
	else if (Optimisation == VTK_OPTIMISATION_SPARSECHOL)
	{
		UseSparseCholeskyOptimisation(dims, OutScalars);
	}
	std::cout << std::endl;
}

std::vector<double> vtkDistanceFieldMRFRegularisation::GetValueTrack()
{
	return ValueTrack;
}

void vtkDistanceFieldMRFRegularisation::ForwardPrediction( vtkDoubleArray* data, vtkDoubleArray* changes, int dim[3] )
{
	int x, y, z;
	int zOffset,yOffset,offset;

	for(z = 0; z < dim[2]; z++)
	{
		zOffset = z*dim[1]*dim[0];
		for(y = 0;y < dim[1]; y++)
		{
			yOffset = y*dim[0] + zOffset;

			for(x = 0;x < dim[0]; x++)
			{
				offset = x + yOffset;
		
				double val    = data->GetValue(offset);
				double change = changes->GetValue(offset);

				double newval = val + PredictionStep * change;

				data->SetValue(offset, newval);
			}
		}
	}
}


void IncreaseValue(vtkIntArray* data, int dims[3], int x, int y, int z, int inc)
{
	int offset = z * dims[0] * dims[1] + y * dims[0] + x;
	if (offset < 0 || offset >= dims[0]*dims[1]*dims[2])
		return;

	int val = data->GetValue(offset) + inc;
	data->SetValue(offset, val);
}

void vtkDistanceFieldMRFRegularisation::UseICMOptimisationWithMemory( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input )
{
	std::cout << "Starting ICM optimisation with memory";
	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		std::cout << " using the double Laplacian prior (Spline)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
	{
		std::cout << " using the voxel difference prior (Membrane)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_UNIFORMGRADIENT)
	{
		std::cout << " using the uniform gradient prior" << std::endl;
	}
	else
	{
		std::cout << " with unknown prior?" << std::endl;
	}


	double spac[3];
	input->GetSpacing(spac);
	double spacing = spac[0];

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

	bool stop = false;

	std::vector<CVoxelID> localVisitList;
	localVisitList.assign(m_VisitOrder.begin(), m_VisitOrder.end());

	for (int it = 0; it < Iterations && !stop; it++)
	{
		// Do not use time to print them all
		if (it < 10 || !(it % 50))
			std::cout << "\rIteration: " << it << " ";

		double sumsq = 0;
		int NChanged = 0;

		int size = localVisitList.size();

		for(int i = 0; i < size; i++)
		{
			CVoxelID idx = m_VisitOrder[i];
			double change = 0;
			if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && !BandedICM)
			{
				change = LocalDoubleLaplacianFiltering(OutScalars, dims, idx);
			}
			else if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && BandedICM)
			{
				change = LocalBandedDoubleLaplacianFiltering(OutScalars, dims, idx);
			}
			else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
			{
				change = LocalSmoothingFiltering(OutScalars, dims, idx);
			}
			else
			{
				change = LocalUniformGradientFiltering(OutScalars, dims, idx, spacing);
			}

			sumsq += change * change;

			// Mark existing band
			int x = idx.id[0];
			int y = idx.id[1];
			int z = idx.id[2];
			IncreaseValue(ICMUpdateVolume, dims, x, y, z, 1000);

			if (change > MinRMS * 0.001)
			{
				NChanged++;
				for (int dx = -1; dx <= 1; dx++)
					for (int dy = -1; dy <= 1; dy++)
						for (int dz = -1; dz <= 1; dz++)
							IncreaseValue(ICMUpdateVolume, dims, x+dx, y+dy, z+dz,1);
				// If you have changed then you are at least 1001
			}
		}
		sumsq /= size;
		double RMS = sqrt(sumsq);

		localVisitList.clear();

		int x, y, z;
		int zOffset,yOffset,offset;

		for(z = 0; z < dims[2]; z++)
		{
			zOffset = z*dims[1]*dims[0];
			for(y = 0;y < dims[1]; y++)
			{
				yOffset = y*dims[0] + zOffset;
				for(x = 0;x < dims[0]; x++)
				{
					offset = x + yOffset;
					if (ICMUpdateVolume->GetValue(offset) > 1000) // only keep the ones in the band and have been updated or a neighbour has been updated
					{
						CVoxelID vid;
						vid.id[0] = x;
						vid.id[1] = y;
						vid.id[2] = z;

						localVisitList.push_back(vid);
						ICMUpdateVolume->SetValue(offset, 0);
					}
				}
			}
		}
		if (it < 10 || !(it % 50))
			std::cout << "RMS (change) " << RMS << " changed: " << NChanged << " !changed: " << size-NChanged << " visits " << localVisitList.size() << "     ";
		if (RMS < MinRMS)
		{
			stop = true;
			std::cout << "Stopping since RMS is less than RMS threshold of " << MinRMS << std::endl;
		}
		if (localVisitList.size() == 0)
		{
			stop = true;
			std::cout << "Stopping since visitlist is empty"  << std::endl;
		}
	}

	std::cout << std::endl;
	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;
}



double vtkDistanceFieldMRFRegularisation::InnerLoopThread( int start, int stop, vtkDoubleArray * OutScalars, int * dims )
{
	double sumsq = 0;
	int NChanged = 0;
	for(int i = start; i < stop; i++)
	{
		CVoxelID idx = m_VisitOrder[i];
		double change = 0;
		if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && !BandedICM)
		{
			change = LocalDoubleLaplacianFiltering(OutScalars, dims, idx);
		}
		else if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && BandedICM)
		{
			change = LocalBandedDoubleLaplacianFiltering(OutScalars, dims, idx);
		}
		else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
		{
			change = LocalSmoothingFiltering(OutScalars, dims, idx);
		}
		//else
		//{
		//	change = LocalUniformGradientFiltering(OutScalars, dims, idx, spacing);
		//}

		sumsq += change * change;
		//if (change > MinRMS)
		//{
		//	NChanged++;
		//}
		//if (change > 1000 || change < -1000)
		//{
		//	std::cerr << "Something wrong with change " << change << std::endl;
		//}

	}
	//sumsq /= size;
	//double RMS = sqrt(sumsq);
	return sumsq;
}


double vtkDistanceFieldMRFRegularisation::InnerLoopWithConvolutionBorderVersionThread( ThreadArgs Targs, bool topbottom )
{
	double sumsq = 0;
	int NChanged = 0;

	// Take entire plane
	if (topbottom)
	{
		for (int z = Targs.start; z < Targs.stop; z++)
		{
			for(int y = 0; y < Targs.dims[1]; y++)
			{
				for(int x = 0; x < Targs.dims[0]; x++)
				{
					double change = 0;
					if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
					{
						change = LocalDoubleLaplacianFilteringUsingConvolutionAndBorderCheck(Targs.OutScalars, Targs.TempOut, Targs.dims, x, y, z);
						sumsq += change * change;
					}

				}
			}
		}
	}
	else // just take the edges
	{
		std::vector<int> ypos;
		std::vector<int> xpos;
		ypos.push_back(0);
		ypos.push_back(1);
		ypos.push_back(Targs.dims[1]-2);
		ypos.push_back(Targs.dims[1]-1);
		xpos.push_back(0);
		xpos.push_back(1);
		xpos.push_back(Targs.dims[0]-2);
		xpos.push_back(Targs.dims[0]-1);
		for (int z = Targs.start; z < Targs.stop; z++)
		{
			for (int yp = 0; yp < ypos.size(); yp++)
			{
				int  y = ypos[yp];
				for (int x = 0; x < Targs.dims[0]; x++)
				{
					double change = 0;
					if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
					{
						change = LocalDoubleLaplacianFilteringUsingConvolutionAndBorderCheck(Targs.OutScalars, Targs.TempOut, Targs.dims, x, y, z);
						sumsq += change * change;
					}

				}
			}
			for (int xp = 0; xp < xpos.size(); xp++)
			{
				int x  = xpos[xp];
				for (int y = 2; y < Targs.dims[1]-2; y++)
				{
					double change = 0;
					if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
					{
						change = LocalDoubleLaplacianFilteringUsingConvolutionAndBorderCheck(Targs.OutScalars, Targs.TempOut, Targs.dims, x, y, z);
						sumsq += change * change;
					}

				}
			}
		}
	}

	return sumsq;
}


// Each loop take one or a few lines in the image
double vtkDistanceFieldMRFRegularisation::InnerLoopWithConvolutionThread( ThreadArgs Targs )
{
	double sumsq = 0;
	int NChanged = 0;

	for (int z = Targs.start; z < Targs.stop; z++)
	{
		for(int y = 2; y < Targs.dims[1]-2; y++)
		{
			for(int x = 2; x < Targs.dims[0]-2; x++)
			{
				double change = 0;
				if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
				{
					change = LocalDoubleLaplacianFilteringUsingConvolution(Targs.OutScalars, Targs.TempOut, Targs.dims, x, y, z);
					sumsq += change * change;
				}

			}
		}
	}
	return sumsq;
}


double vtkDistanceFieldMRFRegularisation::InnerLoopWithSeperateOutputThread( ThreadArgs Targs )
{
	double sumsq = 0;
	int NChanged = 0;
	for(int i = Targs.start; i < Targs.stop; i++)
	{
		CVoxelID idx = m_VisitOrder[i];
		double change = 0;
		if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && !BandedICM)
		{
			change = LocalDoubleLaplacianFilteringWithSeperateOutput(Targs.OutScalars, Targs.TempOut, Targs.dims, idx);
		}
		//else if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && BandedICM)
		//{
		//	change = LocalBandedDoubleLaplacianFiltering(Targs.OutScalars, Targs.dims, idx);
		//}
		//else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
		//{
		//	change = LocalSmoothingFiltering(Targs.OutScalars, Targs.dims, idx);
		//}


		sumsq += change * change;
	}
	return sumsq;
}


void vtkDistanceFieldMRFRegularisation::UseICMOptimisationWithConvolutionMultithreaded( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input )
{
	std::cout << "Starting ICM (convolution) optimisation";
	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		std::cout << " using the double Laplacian prior (Spline) (" << NumThreads << " threads)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
	{
		std::cout << " using the voxel difference prior (Membrane) (" << NumThreads << " threads)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_UNIFORMGRADIENT)
	{
		std::cout << " using the uniform gradient prior (" << NumThreads << " threads)" << std::endl;
	}
	else
	{
		std::cout << " with unknown prior?" << std::endl;
	}

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

	vtkDoubleArray *TempOutput = vtkDoubleArray::New();
	TempOutput->SetNumberOfValues(OutScalars->GetNumberOfTuples());
	for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		TempOutput->SetValue(i, OutScalars->GetValue(i));
	}

	bool stop = false;
	for (int it = 0; it < Iterations && !stop; it++)
	{
		if (it < 10 || !(it % 500))
			std::cout << "\rIteration: " << it << " ";

		const int parts = NumThreads;

		int zmax = dims[2]-4;
		std::vector<int>bnd = CGeneralUtils::bounds(parts, zmax);
		std::vector<std::future<double> > tt;

		//Launch part threads
		// Take care of the main part of the volume
		for (int i = 0; i < parts; ++i) 
		{
			ThreadArgs Targs;
			Targs.start = bnd[i]+2;
			Targs.stop = bnd[i+1]+2;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;

			tt.push_back(std::async(&vtkDistanceFieldMRFRegularisation::InnerLoopWithConvolutionThread, this, Targs));
		}

		int size = dims[0] * dims[1] * dims[2];
		double sumSQ = 0;

		//// Main loop takes the sides
		//{
		//	ThreadArgs Targs;
		//	Targs.start = 2;
		//	Targs.stop = dims[2]-2;
		//	Targs.OutScalars = OutScalars;
		//	Targs.TempOut = TempOutput;
		//	Targs.dims = dims;

		//	sumSQ += InnerLoopWithConvolutionBorderVersionThread(Targs, false);
		//}

		// main loop takes the top
		{
			ThreadArgs Targs;
			Targs.start = 0;
			Targs.stop = 2;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;
			sumSQ += InnerLoopWithConvolutionBorderVersionThread(Targs, true);
		}

		// Gather the remaining threads
		for (int i = 0; i < tt.size(); i++)
		{
			sumSQ += tt[i].get();
		}

		tt.clear();
		// Launch again to take sides
		for (int i = 0; i < parts; ++i) 
		{
			ThreadArgs Targs;
			Targs.start = bnd[i]+2;
			Targs.stop = bnd[i+1]+2;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;

			tt.push_back(std::async(&vtkDistanceFieldMRFRegularisation::InnerLoopWithConvolutionBorderVersionThread, this, Targs, false));
		}
		// main loop takes the bottom
		{
			ThreadArgs Targs;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;
			Targs.start = dims[2]-2;
			Targs.stop = dims[2];
			sumSQ += InnerLoopWithConvolutionBorderVersionThread(Targs, true);
		}

		// Gather the remaining threads again
		for (int i = 0; i < tt.size(); i++)
		{
			sumSQ += tt[i].get();
		}


		sumSQ /= size;

		double RMS = sqrt(sumSQ);

		if (RMS <= 0 || RMS > 10000000)
			std::cout << it << " Something wrong with RMS " << RMS << std::endl;

		if (it < 10 || !(it % 500))
			std::cout << "RMS (change) " << RMS <<  "      ";

		if (RMS < MinRMS)
		{
			stop = true;
			std::cout << "Stopping since RMS is less than RMS threshold of " << MinRMS << std::endl;
		}

		for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
		{
			OutScalars->SetValue(i, TempOutput->GetValue(i));
		}
	}

	TempOutput->Delete();
	std::cout << std::endl;
	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;
}


void vtkDistanceFieldMRFRegularisation::UseICMOptimisationWithConvolution( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input )
{
	std::cout << "Starting ICM optimisation";
	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		std::cout << " using the double Laplacian prior (Spline) (" << NumThreads << " threads)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
	{
		std::cout << " using the voxel difference prior (Membrane) (" << NumThreads << " threads)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_UNIFORMGRADIENT)
	{
		std::cout << " using the uniform gradient prior (" << NumThreads << " threads)" << std::endl;
	}
	else
	{
		std::cout << " with unknown prior?" << std::endl;
	}

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

	vtkDoubleArray *TempOutput = vtkDoubleArray::New();
	TempOutput->SetNumberOfValues(OutScalars->GetNumberOfTuples());
	for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		TempOutput->SetValue(i, OutScalars->GetValue(i));
	}

	bool stop = false;
	for (int it = 0; it < Iterations && !stop; it++)
	{
		//if (it == 0)
		//	DebugVolume.resize(OutScalars->GetNumberOfTuples(), 0);

		if (it < 10 || !(it % 500))
			std::cout << "\rIteration: " << it << " ";

		//const int parts = NumThreads;

		//std::vector<int>bnd = CGeneralUtils::bounds(parts, size);
		//std::vector<std::future<double> > tt;

		////Launch parts-1 threads
		//for (int i = 0; i < parts - 1; ++i) 
		//{
		//	ThreadArgs Targs;
		//	Targs.start = bnd[i];
		//	Targs.stop = bnd[i+1];
		//	Targs.OutScalars = OutScalars;
		//	Targs.TempOut = TempOutput;
		//	Targs.dims = dims;

		//	tt.push_back(std::async(&vtkDistanceFieldMRFRegularisation::InnerLoopWithSeperateOutputThread, this, Targs));
		//}

		//double sumSQ = 0;
		////Use the main thread to do part of the work !!!
		//for (int i = parts - 1; i < parts; ++i) 
		//{
		//	ThreadArgs Targs;
		//	Targs.start = bnd[i];
		//	Targs.stop = bnd[i+1];
		//	Targs.OutScalars = OutScalars;
		//	Targs.TempOut = TempOutput;
		//	Targs.dims = dims;

		//	sumSQ += InnerLoopWithSeperateOutputThread(Targs);
		//}
		//for (int i = 0; i < tt.size(); i++)
		//{
		//	sumSQ += tt[i].get();
		//}


		int size = dims[0] * dims[1] * dims[2];
		double sumSQ = 0;

		// First take the main part (minus the borders)
		{
			ThreadArgs Targs;
			Targs.start = 2;
			Targs.stop = dims[2]-2;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;

			sumSQ += InnerLoopWithConvolutionThread(Targs);
		}
		// now take the sides
		{
			ThreadArgs Targs;
			Targs.start = 2;
			Targs.stop = dims[2]-2;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;

			sumSQ += InnerLoopWithConvolutionBorderVersionThread(Targs, false);
		}
		// now take the top and bottom
		{
			ThreadArgs Targs;
			Targs.start = 0;
			Targs.stop = 2;
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;
			sumSQ += InnerLoopWithConvolutionBorderVersionThread(Targs, true);

			Targs.start = dims[2]-2;
			Targs.stop = dims[2];
			sumSQ += InnerLoopWithConvolutionBorderVersionThread(Targs, true);
		}

		sumSQ /= size;

		double RMS = sqrt(sumSQ);

		if (RMS <= 0 || RMS > 10000000)
			std::cout << it << " Something wrong with RMS " << RMS << std::endl;

		if (it < 10 || !(it % 500))
			std::cout << "RMS (change) " << RMS <<  "      ";

		if (RMS < MinRMS)
		{
			stop = true;
			std::cout << "Stopping since RMS is less than RMS threshold of " << MinRMS << std::endl;
		}

		for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
		{
			OutScalars->SetValue(i, TempOutput->GetValue(i));
		}
		//// Check debug volume
		//if (it == 0)
		//{
		//	std::cout << "Checking debug volume" << std::endl;
		//	for (int i = 0; i < DebugVolume.size(); i++)
		//	{
		//		if (DebugVolume[i] != 1)
		//		{ 
		//			int xt,yt,zt;
		//			CGeneralUtils::GetXYZfromOffset(dims, i, xt, yt, zt);
		//			std::cout << "Debug volume problem at " << i << " (" << xt << ", " << yt << ", " << zt << ") : " << DebugVolume[i] << std::endl;
		//		}
		//	}
		//}
	}

	TempOutput->Delete();
	std::cout << std::endl;
	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;

}


void vtkDistanceFieldMRFRegularisation::UseICMOptimisationWithSeperateOutputMultithreaded( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input )
{
	std::cout << "Starting ICM optimisation";
	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		std::cout << " using the double Laplacian prior (Spline) (" << NumThreads << " threads)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
	{
		std::cout << " using the voxel difference prior (Membrane) (" << NumThreads << " threads)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_UNIFORMGRADIENT)
	{
		std::cout << " using the uniform gradient prior (" << NumThreads << " threads)" << std::endl;
	}
	else
	{
		std::cout << " with unknown prior?" << std::endl;
	}

	int size = m_VisitOrder.size();

	//int TRx = dims[0] / 2;
	//int TRy = dims[0] / 2;
	//int TRz = dims[0] / 2;

	double spac[3];
	input->GetSpacing(spac);
	double spacing = spac[0];

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

	vtkDoubleArray *TempOutput = vtkDoubleArray::New();
	TempOutput->SetNumberOfValues(OutScalars->GetNumberOfTuples());
	for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		TempOutput->SetValue(i, OutScalars->GetValue(i));
	}
//	DebugVolume.resize(OutScalars->GetNumberOfTuples(), 0);


	bool stop = false;
	for (int it = 0; it < Iterations && !stop; it++)
	{
		if (it < 10 || !(it % 500))
			std::cout << "\rIteration: " << it << " ";

		const int parts = NumThreads;

		std::vector<int>bnd = CGeneralUtils::bounds(parts, size);
		std::vector<std::future<double> > tt;

		//Launch parts-1 threads
		for (int i = 0; i < parts - 1; ++i) 
		{
			ThreadArgs Targs;
			Targs.start = bnd[i];
			Targs.stop = bnd[i+1];
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;

			tt.push_back(std::async(&vtkDistanceFieldMRFRegularisation::InnerLoopWithSeperateOutputThread, this, Targs));
		}

		double sumSQ = 0;
		//Use the main thread to do part of the work !!!
		for (int i = parts - 1; i < parts; ++i) 
		{
			ThreadArgs Targs;
			Targs.start = bnd[i];
			Targs.stop = bnd[i+1];
			Targs.OutScalars = OutScalars;
			Targs.TempOut = TempOutput;
			Targs.dims = dims;

			sumSQ += InnerLoopWithSeperateOutputThread(Targs);
		}
		for (int i = 0; i < tt.size(); i++)
		{
			sumSQ += tt[i].get();
		}


		//double sumSQ = 0;
		//{
		//	ThreadArgs Targs;
		//	Targs.start = 0;
		//	Targs.stop = size;
		//	Targs.OutScalars = OutScalars;
		//	Targs.TempOut = TempOutput;
		//	Targs.dims = dims;

		//	sumSQ += InnerLoopWithSeperateOutputThread(Targs);
		//}
		//

		sumSQ /= size;
		double RMS = sqrt(sumSQ);

		if (RMS <= 0 || RMS > 10000000)
			std::cout << it << " Something wrong with RMS " << RMS << std::endl;

		if (it < 10 || !(it % 500))
			std::cout << "RMS (change) " << RMS <<  "      ";

		if (RMS < MinRMS)
		{
			stop = true;
			std::cout << "Stopping since RMS is less than RMS threshold of " << MinRMS << std::endl;
		}

		//double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
		//double StartESpline2 = EnergyComputer.TotalSplineEnergy(TempOutput);
		//std::cout << StartESpline << " " << StartESpline2 << std::endl;
		
		//		OutScalars->DeepCopy(TempOutput);
		for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
		{
			OutScalars->SetValue(i, TempOutput->GetValue(i));
		}

		//StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
		//StartESpline2 = EnergyComputer.TotalSplineEnergy(TempOutput);
		//std::cout << StartESpline << " " << StartESpline2 << std::endl;
	}

	TempOutput->Delete();
	std::cout << std::endl;
	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;

	//// Check debug volume
	//for (int i = 0; i < DebugVolume.size(); i++)
	//{
	//	if (DebugVolume[i] != 1)
	//	{
	//		std::cout << "Debug volume problem at " << i << std::endl;
	//	}
	//}
}

void vtkDistanceFieldMRFRegularisation::UseICMOptimisationMultithreaded( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input )
{
	std::cout << "Starting ICM optimisation";
	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		std::cout << " using the double Laplacian prior (Spline)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
	{
		std::cout << " using the voxel difference prior (Membrane)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_UNIFORMGRADIENT)
	{
		std::cout << " using the uniform gradient prior" << std::endl;
	}
	else
	{
		std::cout << " with unknown prior?" << std::endl;
	}

	int size = m_VisitOrder.size();

	int TRx = dims[0] / 2;
	int TRy = dims[0] / 2;
	int TRz = dims[0] / 2;

	double spac[3];
	input->GetSpacing(spac);
	double spacing = spac[0];

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

	bool stop = false;
	for (int it = 0; it < Iterations && !stop; it++)
	{
		if (it < 10 || !(it % 500))
			std::cout << "\rIteration: " << it << " ";

		const int parts = NumThreads;

		std::vector<int>bnd = CGeneralUtils::bounds(parts, size);
		std::vector<std::future<double> > tt;

		//Launch parts-1 threads
		for (int i = 0; i < parts - 1; ++i) 
		{
			tt.push_back(std::async(&vtkDistanceFieldMRFRegularisation::InnerLoopThread, this, bnd[i], bnd[i+1], OutScalars, dims));
		}

		double sumSQ = 0;
		//Use the main thread to do part of the work !!!
		for (int i = parts - 1; i < parts; ++i) 
		{
			sumSQ += InnerLoopThread(bnd[i], bnd[i+1], OutScalars, dims);
		}
		for (int i = 0; i < tt.size(); i++)
		{
			sumSQ += tt[i].get();
		}
		sumSQ /= size;
		double RMS = sqrt(sumSQ);

		if (it < 10 || !(it % 500))
			std::cout << "RMS (change) " << RMS <<  "      ";

		if (RMS < MinRMS)
		{
			stop = true;
			std::cout << "Stopping since RMS is less than RMS threshold of " << MinRMS << std::endl;
		}
	}

	std::cout << std::endl;
	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;
}

void vtkDistanceFieldMRFRegularisation::UseICMOptimisation( int * dims, vtkDoubleArray * OutScalars, vtkImageData* input )
{
	std::cout << "Starting ICM optimisation";
	if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE)
	{
		std::cout << " using the double Laplacian prior (Spline)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
	{
		std::cout << " using the voxel difference prior (Membrane)" << std::endl;
	}
	else if (PriorType == VTK_PRIOR_ENERGY_UNIFORMGRADIENT)
	{
		std::cout << " using the uniform gradient prior" << std::endl;
	}
	else
	{
		std::cout << " with unknown prior?" << std::endl;
	}

	//ChangeVolume = vtkImageData::New();
	//ChangeVolume->CopyStructure(input);
	//ChangeVolume->AllocateScalars();
	//vtkDoubleArray *ChangeScalars = 
	//	vtkDoubleArray::SafeDownCast(ChangeVolume->GetPointData()->GetScalars());

	//if (!ChangeScalars)
	//{
	//	std::cerr << "Something wrong with change volume" << std::endl;
	//	return;
	//}

//	int size = dims[0]*dims[1]*dims[2];
	int size = m_VisitOrder.size();

	//bool DoForwardPrediction = false;

	int TRx = dims[0] / 2;
	int TRy = dims[0] / 2;
	int TRz = dims[0] / 2;

	double spac[3];
	input->GetSpacing(spac);
	double spacing = spac[0];

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

	//EnergyComputer.Test(OutScalars, 0, 0, 0);
	//EnergyComputer.Test(OutScalars, 10, 15, 20);
	//EnergyComputer.Test(OutScalars, 10, 0, 0);
	//EnergyComputer.Test(OutScalars, 0, 15, 0);
	//EnergyComputer.Test(OutScalars, 0, 0, 20);

	//EnergyComputer.Test(OutScalars, 0);
	//EnergyComputer.Test(OutScalars, 10000);
	//EnergyComputer.Test(OutScalars, 20000);
	//EnergyComputer.Test(OutScalars, 30000);

//	int lastForward = 10;
	bool stop = false;
//	double lowestRMS = VTK_FLOAT_MAX;
	for (int it = 0; it < Iterations && !stop; it++)
	{
		double trackval = GetValue(OutScalars, dims, TRx, TRy, TRz);
		ValueTrack.push_back(trackval);

		// Do not use time to print them all
		if (it < 10 || !(it % 50))
			std::cout << "\rIteration: " << it << " ";

		double sumsq = 0;
		int NChanged = 0;
		for(int i = 0; i < size; i++)
		{
			//		int idx = m_VisitOrder[i];
			//		OutScalars->SetValue(idx, InScalars->GetValue(idx));

			CVoxelID idx = m_VisitOrder[i];
			double change = 0;
			if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && !BandedICM)
			{
				change = LocalDoubleLaplacianFiltering(OutScalars, dims, idx);
			}
			else if (PriorType == VTK_PRIOR_ENERGY_DIFF_LAPLACE && BandedICM)
			{
				change = LocalBandedDoubleLaplacianFiltering(OutScalars, dims, idx);
			}
			else if (PriorType == VTK_PRIOR_ENERGY_DIFF_VOXEL)
			{
				change = LocalSmoothingFiltering(OutScalars, dims, idx);
			}
			else
			{
				change = LocalUniformGradientFiltering(OutScalars, dims, idx, spacing);
			}

//			SetValue(ChangeScalars, dims, idx.id[0], idx.id[1], idx.id[2], change);
			sumsq += change * change;

			if (change > MinRMS)
			{
				NChanged++;
			}
		}
		sumsq /= size;
		double RMS = sqrt(sumsq);

		// Do not use time to print them all
		if (it < 10 || !(it % 50))
			std::cout << "RMS (change) " << RMS << " changed: " << NChanged << " !changed: " << size-NChanged << "      ";
	
		if (RMS < MinRMS)
		{
			stop = true;
			std::cout << "Stopping since RMS is less than RMS threshold of " << MinRMS << std::endl;
		}
		//if (DoForwardPrediction && RMS < 0.03 && lastForward < 0 && (Iterations-it) > 10 && RMS < lowestRMS)
		//{
		//	std::cout << "Forward prediction with " << PredictionStep << " steps" << std::endl;
		//	ForwardPrediction(OutScalars, ChangeScalars, dims);
		//	lastForward	 = 11;
		//	PredictionStep *= 2.0;
		//	if (PredictionStep > 60)
		//		PredictionStep = 60;
		//}
		//lowestRMS = std::min(lowestRMS, RMS);
		//lastForward--;
	}

	std::cout << std::endl;
	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;
}

void vtkDistanceFieldMRFRegularisation::UseConjugateGradientOptimisation( int * dims, vtkDoubleArray * OutScalars )
{
	std::cout << "Starting Conjugate Gradient optimisation" << std::endl;

	CMRFEnergyComputer EnergyComputer(dims, OrgScalars, LocalWeights, GlobalBeta);
	double StartEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	double StartESpline = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl << "Starting energy membrane " << StartEMemb << " spline " << StartESpline << std::endl << std::endl;

//	EnergyComputer.Test(OutScalars, 0, 0, 0);
//	EnergyComputer.Test(OutScalars, 10, 15, 20);
//	EnergyComputer.Test(OutScalars, 10, 0, 0);
//	EnergyComputer.Test(OutScalars, 0, 15, 0);
//	EnergyComputer.Test(OutScalars, 0, 0, 20);

//	EnergyComputer.TestGradAndTot(OutScalars);

	CMRFCostFunction costFunction(dims, OutScalars->GetNumberOfTuples(), OrgScalars, LocalWeights, GlobalBeta, PriorType);

	vnl_vector<double> d(OutScalars->GetNumberOfTuples());
	for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		d[i] = OutScalars->GetValue(i);
	}

	vnl_conjugate_gradient cg(costFunction); 
	cg.set_f_tolerance(CGTolerance);
	cg.set_g_tolerance(CGTolerance);
	cg.set_x_tolerance(1e-5);
	cg.set_max_function_evals(1000);
//	cg.set_max_function_evals(10);
	cg.minimize(d);
	std::cout << std::endl;
	cg.diagnose_outcome();

	//vnl_lbfgs lbfgs(costFunction);
	//lbfgs.set_max_function_evals(500);
	//lbfgs.minimize(d);

	for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		OutScalars->SetValue(i, d[i]);
	}

	double FinalEMemb = EnergyComputer.TotalMembraneEnergy(OutScalars);
	FinalSplineEnergy = EnergyComputer.TotalSplineEnergy(OutScalars);
	std::cout << std::endl <<  "Final energy membrane " << FinalEMemb << " spline " << FinalSplineEnergy << std::endl;
}



void vtkDistanceFieldMRFRegularisation::CreateLookup(int * dims, std::vector<int> &Lookup, std::vector<vgl_point_3d<int> >& PosLU)
{
	if (Lookup.size() > 0)
	{
		Lookup.clear();
		PosLU.clear();
	}

	// -x
	Lookup.push_back(-1);
	// +x
	Lookup.push_back( 1);

	// -y
	Lookup.push_back(-dims[0]);

	// y
	Lookup.push_back( dims[0]);

	// -z
	Lookup.push_back(-dims[1] * dims[0]);

	// z
	Lookup.push_back( dims[1] * dims[0]);

	PosLU.push_back(vgl_point_3d<int>(-1, 0, 0));
	PosLU.push_back(vgl_point_3d<int>( 1, 0, 0));
	PosLU.push_back(vgl_point_3d<int>( 0,-1, 0));
	PosLU.push_back(vgl_point_3d<int>( 0, 1, 0));
	PosLU.push_back(vgl_point_3d<int>( 0, 0, -1));
	PosLU.push_back(vgl_point_3d<int>( 0, 0,  1));
}

void vtkDistanceFieldMRFRegularisation::FindMatrixEntriesForNeighbours( int * dims, std::vector<int> &Lookup, std::vector<vgl_point_3d<int> >& PosLU, 
																	   int x, int y, int z, int offset, std::vector<CMatrixEntry>& entries, std::vector<double> &b)
{
	// alpha = 1: only original distance. Very high local sampling density
	// alpha = 0: pure smoothing. Very low local sampling density
	double localAlpha = GetValue(LocalWeights, dims, x, y, z) * GlobalBeta;
	double orgvalue   = GetValue(OrgScalars, dims, x, y, z);

	// Trivial solution to local equation
	// di = di_org
	if (localAlpha >= 1.0)
	{
		CMatrixEntry entry;
		entry.i = offset;
		entry.j = offset;
		entry.v = 1;
		entries.push_back(entry);

		b[offset] = orgvalue;
		return;
	}

	int N = 0;
	for (unsigned int i = 0; i < PosLU.size(); i++)
	{
		vgl_point_3d<int> p(PosLU[i].x() + x, PosLU[i].y() + y, PosLU[i].z() + z);

		if (p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < dims[0] && p.y() < dims[1] && p.z() < dims[2])
		{
			N++;

			int tid = offset + Lookup[i];
			
			CMatrixEntry entry;
			entry.i = offset;
			entry.j = tid;
			entry.v = -1;
			
			// Only lower triangular
			if (entry.i >= entry.j)
				entries.push_back(entry);
		}
	}
	if (N == 0)
	{
		std::cerr << "N == 0 in FindMatrixEntriesForNeighbours";
		return;
	}

	CMatrixEntry entry;
	entry.i = offset;
	entry.j = offset;
	entry.v = N / (1-localAlpha);
	entries.push_back(entry);

	b[offset] = N * localAlpha / (1-localAlpha) * orgvalue;
}

void vtkDistanceFieldMRFRegularisation::CheckMatrixSymmetry( std::vector<CMatrixEntry>& matrix, std::vector<double>& b )
{
	//! very very slow
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		CMatrixEntry entry = matrix[i];
		
		// Check diag
		if (entry.i != entry.j)
		{
			//! Now find transpose
			bool found = false;
			for (unsigned int j = 0; j < matrix.size() && !found; j++)
			{
				if (entry.i == matrix[j].j && entry.j == matrix[j].i)
				{
					if (entry.v == matrix[j].v)
					{
						found = true;
					}
					else
					{
						std::cout << "values wrong at " << entry.i << " " << entry.j << std::endl;
					}
				}
			}
			if (!found)
			{
				std::cout << "No transpose at " << entry.i << " " << entry.j << std::endl;
			}
		}
	}
}

void vtkDistanceFieldMRFRegularisation::ReadResult(vtkDoubleArray * OutScalars)
{
	std::string iname = "C:\\rrplocal\\data\\IMM\\MRFSurf\\x.mtx";
	std::ifstream fist(iname.c_str());

	std::string t;
	std::getline(fist, t);
	std::getline(fist, t);

	int N = OutScalars->GetNumberOfTuples();
	for (int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		double v = 0;
		fist >> v;
		if (fist.fail())
		{
			std::cerr << "read error" << std::endl;
		}
		OutScalars->SetValue(i, v);
	}
}

void vtkDistanceFieldMRFRegularisation::WriteMatrixAndB(std::vector<CMatrixEntry>& matrix, std::vector<double>& b)
{
	std::string oname1 = "C:\\rrplocal\\data\\IMM\\MRFSurf\\matrix.csv";
	std::string oname2 = "C:\\rrplocal\\data\\IMM\\MRFSurf\\b.txt";

	std::ofstream ost1(oname1.c_str());
	
//	ost1 << matrix.size() << std::endl;
	ost1 << std::setprecision(8);
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		if (matrix[i].i >= matrix[i].j)
		{
			ost1 << matrix[i].i << ", " << matrix[i].j << ", " << matrix[i].v << std::endl;
		}
	}

	std::ofstream ost2(oname2.c_str());
	ost2 << b.size() << std::endl;

	for (unsigned int i = 0; i < b.size(); i++)
	{
		ost2 << b[i] << std::endl;
	}
}

void vtkDistanceFieldMRFRegularisation::CreateSparseMatrix( int * dims, std::vector<CMatrixEntry> &matrix, std::vector<double>& b )
{
	//! Lookup for neighbour voxels using offsets
	std::vector<int> Lookup;

	//! Lookup using x, y,z 
	std::vector<vgl_point_3d<int> > PosLU;

	CreateLookup(dims, Lookup, PosLU);

	int x, y, z;
	int zOffset,yOffset,offset;

	for(z = 0; z < dims[2]; z++)
	{
		zOffset = z * dims[1] * dims[0];
		for(y = 0; y < dims[1]; y++)
		{
			yOffset = y * dims[0] + zOffset;

			for(x = 0; x < dims[0]; x++)
			{
				offset = x + yOffset;

				FindMatrixEntriesForNeighbours(dims, Lookup, PosLU, x, y, z, offset, matrix, b);
			}
		}
	}
	std::cout << "Sparse matrix with " << matrix.size() << " elements" << std::endl;
}

void vtkDistanceFieldMRFRegularisation::UseSparseCholeskyOptimisation( int * dims, vtkDoubleArray * OutScalars )
{
	#ifdef USECHOLMOD
	std::string oname = "C:\\rrplocal\\data\\IMM\\MRFSurf\\testmatrix.mtx";

	int size = dims[0] * dims[1] * dims[2];

	std::vector<CMatrixEntry> matrix;
	std::vector<double> bside(size);

	std::cout << "Creating matrices" << std::endl;
	CreateSparseMatrix(dims, matrix, bside);

	cholmod_common c ;
	cholmod_start (&c) ;		
	    
	int nrow   = OutScalars->GetNumberOfTuples();
	int ncol   = OutScalars->GetNumberOfTuples();
	int nnz    = matrix.size(); // non zero elements
	int stype = 1; // Lower triangular
	int xtype = 1; // Real
	cholmod_triplet *ATrip = cholmod_allocate_triplet(nrow, ncol, nnz, stype, xtype, &c);

	int *is = (int*)ATrip->i;
	int *js = (int*)ATrip->j;
	double *xs = (double*)ATrip->x;

	if (is == NULL || js == NULL || xs == NULL)
	{
		std::cerr << "Could not allocate triplet matrix" << std::endl;
		return;
	}

	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		is[i] = matrix[i].i;
		js[i] = matrix[i].j;
		xs[i] = matrix[i].v;
	}
	ATrip->nnz = nnz;

	int isok = cholmod_check_triplet(ATrip, &c);


	cholmod_sparse *A = cholmod_triplet_to_sparse(ATrip, 0, &c);
	cholmod_free_triplet(&ATrip, &c);

	//FILE *f3 = fopen(oname.c_str(), "w");
	//cholmod_write_sparse(f3, A, NULL, NULL, &c);
	//fclose(f3);


	cholmod_dense *b  = cholmod_allocate_dense(nrow, 1, nrow, xtype, &c);
	double *xs1 = (double*)b->x;
	for (unsigned int i = 0; i < bside.size(); i++)
	{
		xs1[i] = bside[i];
	}
	
	cholmod_print_sparse (A, "A", &c) ;		    /* print the matrix */
	//if (A == NULL || A->stype == 0)		    /* A must be symmetric */
	//{
	//	std::cout << "Something wrong with A" << std::endl;
	//	cholmod_free_sparse (&A, &c) ;
	//	cholmod_finish (&c) ;
	//	return (0) ;
	//}
	cholmod_print_dense(b, "b", &c);

	std::cout << "Factorising" << std::endl;
	cholmod_factor *L = cholmod_analyze (A, &c);
	cholmod_factorize (A, L, &c) ;		    

	std::cout << "Solving" << std::endl;
	cholmod_dense *x = cholmod_solve (CHOLMOD_A, L, b, &c);

	cholmod_print_dense(x, "x", &c);

	double *xs2 = (double*)x->x;

	for ( int i = 0; i < OutScalars->GetNumberOfTuples(); i++)
	{
		OutScalars->SetValue(i, xs2[i]);
	}

	//FILE *f3 = fopen(oname.c_str(), "w");
	//cholmod_write_dense(f3, x, NULL, &c);


	////	cholmod_dense *r  = cholmod_copy_dense (b, &c) ;		    /* r = b */
	////cholmod_sdmult (A, 0, m1, one, x, r, &c) ;	    /* r = r-Ax */
	////printf ("norm(b-Ax) %8.1e\n",
	////	cholmod_norm_dense (r, 0, &c)) ;	    /* print norm(r) */
	
	cholmod_free_factor (&L, &c) ;		    /* free matrices */
	cholmod_free_sparse (&A, &c) ;

	////	cholmod_free_dense (&r, &c) ;
	cholmod_free_dense (&x, &c) ;
	cholmod_free_dense (&b, &c) ;
	cholmod_finish (&c) ;			    /* finish CHOLMOD */

	//fclose(f1);
	//fclose(f2);
	//fclose(f3);


	//CheckMatrixSymmetry(matrix, b);
	//WriteMatrixAndB(matrix, b);
	//ReadResult(OutScalars);
	#else
		std::cerr << "CHOLMOD not included in this build" << std::endl;
	#endif
}

void vtkDistanceFieldMRFRegularisation::SetWeightVol( vtkImageData* WV )
{
	WeightVol = WV;
	LocalWeights = vtkDoubleArray::SafeDownCast(WeightVol->GetPointData()->GetScalars());
	DeleteWeightVol = false;
	UseLocalWeighting = true;
}







CMRFCostFunction::CMRFCostFunction()
{
}

CMRFCostFunction::CMRFCostFunction( int dims[3], int nunknowns,  vtkDoubleArray *OS, vtkDoubleArray *LW, double GB, int EnergyType) : vnl_cost_function(nunknowns)
{
	m_EnergyComputer.SetDimensions(dims);
	m_EnergyComputer.LocalWeights = LW;
	m_EnergyComputer.OrgScalars = OS;
	m_EnergyComputer.GlobalBeta = GB;
	m_EnergyType = EnergyType;
}

double CMRFCostFunction::f( vnl_vector< double > const& x )
{
	double E = 0;
		
	if (m_EnergyType == 1)
		E = m_EnergyComputer.TotalMembraneEnergy(x);
	else
		E = m_EnergyComputer.TotalSplineEnergy(x);

	return E;
}

void CMRFCostFunction::gradf( vnl_vector< double > const &x, vnl_vector< double > &gradient )
{
	//std::cout << "Gradient" << std::endl;
	//fdgradf(x, gradient);
	//std::cout << "Grad end" << std::endl;
	if (m_EnergyType == 1)
		m_EnergyComputer.MembraneEnergyGradient(x, gradient);
	else
		m_EnergyComputer.SplineEnergyGradient(x, gradient);
}

