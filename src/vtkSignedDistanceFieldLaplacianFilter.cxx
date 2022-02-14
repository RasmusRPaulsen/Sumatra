#include "vtkSignedDistanceFieldLaplacianFilter.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include <vtkDoubleArray.h>
#include <vtkpointdata.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>

//vtkCxxRevisionMacro(vtkSignedDistanceFieldLaplacianFilter, "$Revision: 1.11 $");
vtkStandardNewMacro(vtkSignedDistanceFieldLaplacianFilter);

vtkSignedDistanceFieldLaplacianFilter::vtkSignedDistanceFieldLaplacianFilter()
{
	Iterations = 1;
}

void vtkSignedDistanceFieldLaplacianFilter::VisualiseCoefficients()
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

void vtkSignedDistanceFieldLaplacianFilter::CreateCoefficients()
{
	m_FilterCoefficients.clear();
	
	double scale = -m_CoeffCube[2][2][2];
	if (scale == 0)
	{
		std::cerr << "Something is wrong" << std::endl;
		return;
	}

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


void vtkSignedDistanceFieldLaplacianFilter::LocalLaplaceCoefficients(int x, int y, int z, double fac)
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

void vtkSignedDistanceFieldLaplacianFilter::CreateLaplaceCoefficients()
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

void vtkSignedDistanceFieldLaplacianFilter::CreateVisitOrder(int dim[3])
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

double GetValue(vtkDoubleArray* data, int dims[3], int x, int y, int z)
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

void SetValue(vtkDoubleArray* data, int dims[3], int x, int y, int z, double val)
{
	int offset = z * dims[0] * dims[1] + y * dims[0] + x;
	data->SetValue(offset, val);
}

void vtkSignedDistanceFieldLaplacianFilter::LocalDoubleLaplacianFiltering(vtkDoubleArray* data, int dims[3], CVoxelID idx)
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

	SetValue(data, dims, x, y, z, val);
}


void vtkSignedDistanceFieldLaplacianFilter::LocalSmoothingFiltering(vtkDoubleArray* data,
														   int dims[3],
														   vtkSignedDistanceFieldLaplacianFilter::CVoxelID idx)
{
	int x = idx.id[0];
	int y = idx.id[1];
	int z = idx.id[2];

	//// Handle borders - badly
	//if (x == 0 || x == dims[0]-1 || y == 0 || y == dims[1]-1 || z == 0 || z == dims[2]-1)
	//{
	//	return;
	//}
	
	double val = GetValue(data, dims, x-1, y  , z   ) +
				 GetValue(data, dims, x+1, y  , z   ) +
				 GetValue(data, dims, x  , y-1, z   ) +
				 GetValue(data, dims, x  , y+1, z   ) +
				 GetValue(data, dims, x  , y  , z-1 ) +
				 GetValue(data, dims, x  , y  , z+1 );

	val /= 6.0;
	
	SetValue(data, dims, x, y, z, val);
}

void vtkSignedDistanceFieldLaplacianFilter::SimpleExecute(vtkImageData* input,
                                                vtkImageData* output)
{
	int i;

	vtkDoubleArray *InScalars = 
		vtkDoubleArray::SafeDownCast(input->GetPointData()->GetScalars());

	vtkDoubleArray *OutScalars = 
		vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

	if (!InScalars || !OutScalars)
	{
		std::cerr << "Something wrong with input or output" << std::endl;
		return;
	}

	int dims[3];
	input->GetDimensions(dims);

	int size = dims[0]*dims[1]*dims[2];
	CreateVisitOrder(dims);

	CreateLaplaceCoefficients();

	// Copy data
	for(i = 0; i < size; i++)
	{
		OutScalars->SetValue(i, InScalars->GetValue(i));
	}

	for (int it = 0; it < Iterations; it++)
	{
		std::cout << "Iteration: " << it << std::endl;
		for(i = 0; i < size; i++)
		{
			//		int idx = m_VisitOrder[i];
			//		OutScalars->SetValue(idx, InScalars->GetValue(idx));

			CVoxelID idx = m_VisitOrder[i];
//			LocalSmoothingFiltering(OutScalars, dims, idx);
			LocalDoubleLaplacianFiltering(OutScalars, dims, idx);
		}
	}
}
