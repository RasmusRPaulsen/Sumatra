#include "vtkExtMisc.h"
#include "vtkExtMisc.h"
#include "GeneralUtils.h"
#include <fstream>
#include <deque>
#include <vector>
#include <iostream>
#include <string>
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCutter.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkPlane.h"
#include "vtkPlaneSource.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPrincipalAxisTransform.h"
#include "vtkSphereSource.h"
#include "vtkSTLReader.h"
#include "vtkOBJReader.h"
#include <vtkOBJExporter.h>
#include <vtkSTLWriter.h>
#include "vtkStripper.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include <vtkTriangle.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkFloatArray.h>
#include <vtkFeatureEdges.h>
#include <sstream>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
//#include "vtkAranzReader.h"
#include "vtkPLYReader.h"
//#include "vtk3DMDTxtReader.h"
#include "vtkPolyDataTextReader.h"
#include "vtkRawPolyWriter.h"
//#include "vtkARANZWriter.h"
#include "vtkPLYWriter.h"
#include "LandmarkCollection.h"
#include "vtkAppendPolyData.h"
#include "vtkXMLImageDataWriter.h"
#include <vtkDenseArray.h>
#include <vtkArrayWriter.h>
#include <vtkArrayReader.h>
#include "vtkOBJWriter.h"
#include <vtkArray.h>
#include <math.h>
#include "vtkLine.h"
#include <vtkCellData.h>


void vtkExtMisc::PrintPolyInfo(vtkPolyData* pd)
{
	std::cout << "Number of points: " << pd->GetNumberOfPoints() << std::endl;
	std::cout << "Number of polys: " << pd->GetNumberOfPolys() << std::endl;
	std::cout << "Number of lines: " << pd->GetNumberOfLines() << std::endl;
	std::cout << "Number of verts: " << pd->GetNumberOfVerts() << std::endl;
};


bool vtkExtMisc::WritePDVTK(vtkPolyData *pd, const std::string& fname, bool Binary)
{
	vtkPolyDataWriter *Writer = vtkPolyDataWriter::New();
 	 Writer->SetInputData(pd);
	 Writer->SetFileName(fname.c_str());
	 if (Binary)
		 Writer->SetFileTypeToBinary();
	 int result = Writer->Write();

	 Writer->Delete();

	return (result == 1);
}

bool vtkExtMisc::WritePDFPointsAsSpheres(vtkPoints * points, double radius, const std::string & fname)
{
	double size = radius;

	if (size == 0)
	{
		// Get an estimate of the size of the pointcloud
		double bounds[6];
		points->GetBounds(bounds);

		double l = sqrt((bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
			(bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
			(bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));

		// Sphere size is a % of bounding box diagonal
		size = l * 0.005;
	}

	vtkAppendPolyData *append = vtkAppendPolyData::New();

	for (int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		vtkSphereSource *sphere = vtkSphereSource::New();
		sphere->SetCenter(points->GetPoint(i));
		sphere->SetRadius(size);
		sphere->SetThetaResolution(20);
		sphere->SetPhiResolution(20);
		sphere->Update();
		append->AddInputData(sphere->GetOutput());

		sphere->Delete();
	}
	append->Update();
	vtkExtMisc::WritePDVTK(append->GetOutput(), fname);

	append->Delete();

	return true;
}

bool vtkExtMisc::WritePDSTL(vtkPolyData *pd, const std::string& fname, bool binary)
{
	vtkSTLWriter *Writer = vtkSTLWriter::New();
 	 Writer->SetInputData(pd);
	 Writer->SetFileName(fname.c_str());
	 if (binary)
		 Writer->SetFileTypeToBinary();

	int result = Writer->Write();

	Writer->Delete();

	return (result == 1);
}

bool vtkExtMisc::WritePDSTLWithNormals(vtkPolyData *pd, const std::string& fname, bool binary)
{
	vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
	norms->SetInputData(pd);
	norms->SplittingOff();
	norms->ConsistencyOn();
	norms->Update();

	vtkSTLWriter *Writer = vtkSTLWriter::New();
	Writer->SetInputData(norms->GetOutput());
	Writer->SetFileName(fname.c_str());
	if (binary)
		Writer->SetFileTypeToBinary();

	int result = Writer->Write();

	Writer->Delete();
	norms->Delete();

	return (result == 1);
}

bool vtkExtMisc::WritePDOBJ(vtkRenderWindow *ren, const std::string& fname)
{
	vtkOBJExporter *Writer = vtkOBJExporter::New();
// 	 Writer->SetReferenceCount(ren);
	 Writer->SetFilePrefix(fname.c_str());
	 Writer->Write();

	 Writer->Delete();

	return true;
}


void vtkExtMisc::WritePDWithNormals(vtkPolyData *pd, const std::string& fname)
{
	vtkPolyDataNormals *norms = vtkPolyDataNormals::New();
	 norms->SetInputData(pd);
	 norms->SplittingOff();
	 norms->ConsistencyOn();
	 norms->Update();

	vtkPolyDataWriter *Writer = vtkPolyDataWriter::New();
 	 Writer->SetInputConnection(norms->GetOutputPort());
	 Writer->SetFileName(fname.c_str());
	 Writer->Write();

	 Writer->Delete();
	 norms->Delete();
}


//! Write a vtkPoints as a polydata with vertex structures
void vtkExtMisc::WritePointsAsVertices(vtkPoints *pts, const std::string& fname)
{
	vtkPolyData *pd = vtkPolyData::New();
	vtkCellArray *verts  = vtkCellArray::New();

	for (int i = 0; i < pts->GetNumberOfPoints(); i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
	}
	pd->SetPoints(pts);
	pd->SetVerts(verts);
	verts->Delete();

	WritePDVTK(pd, fname);

	pd->Delete();
}


//------------------------------------------------------------------------
double** vtkExtMisc::NewMatrix(int rows, int cols) 
{
	double *matrix = new double[rows*cols];
	double **m = new double *[rows];
	for(int i = 0; i < rows; i++) 
	{
		m[i] = &matrix[i*cols];
	}
	return m;
}

//------------------------------------------------------------------------
void vtkExtMisc::DeleteMatrix(double **m) 
{
	delete [] *m;
	delete [] m;
}

//------------------------------------------------------------------------
void vtkExtMisc::ZeroMatrix(double **m, int rows, int cols) 
{
	for(int i = 0; i < rows; i++) 
	{
		for(int j = 0; j < cols; j++) 
		{
			m[i][j] = 0.0;
		}
	}
}

//------------------------------------------------------------------------
void vtkExtMisc::MatrixMultiply(double **a, double **b, double **c,
								int arows, int acols, 
								int brows, int bcols) 
{
	if(acols != brows) 
	{
		return; // acols must equal br otherwise we can't proceed
	}
	
	// c must have size arows*bcols (we assume this)
	
	for(int i = 0; i < arows; i++) 
	{
		for(int j = 0; j < bcols; j++) 
		{
			c[i][j] = 0.0;
			for(int k = 0; k < acols; k++)
			{
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
}

//------------------------------------------------------------------------
// Here it is assumed that a rows >> a cols
// Output matrix is [a rows X a rows]
void vtkExtMisc::LargeCovarianceMatrix(double **a, double **c,
									 int arows, int acols) 
{	
	// c must have size arows*arows (we assume this)
	for(int i = 0; i < arows; i++) 
	{
		for(int j = 0; j < arows; j++) 
		{
			c[i][j] = 0.0;
			for(int k = 0; k < acols; k++)
			{
				c[i][j] += a[i][k]*a[j][k];
			}
		}
	}
}	

//------------------------------------------------------------------------
// Here it is assumed that a rows >> a cols
// Output matrix is [a cols X a cols]
void vtkExtMisc::SmallCovarianceMatrix(double **a, double **c,
									 int arows, int acols) 
{	
	const int s = acols;

	// c must have size acols*acols (we assume this)
	for(int i = 0; i < acols; i++) 
	{
		for(int j = 0; j < acols; j++) 
		{
			c[i][j] = 0.0;
			for(int k = 0; k < arows; k++)
			{
				c[i][j] += a[k][i]*a[k][j];
			}
			c[i][j] /= (s-1);
		}
	}
}	


//------------------------------------------------------------------------
void vtkExtMisc::MatrixTranspose(double **a, double **b, int rows, int cols)
{
	for(int i = 0; i < rows; i++) 
	{
		for(int j = 0; j < cols; j++) 
		{
			double tmp = a[i][j];
			b[i][j] = a[j][i];
			b[j][i] = tmp;
		}
	}
}

//------------------------------------------------------------------------
double* vtkExtMisc::NewVector(int length)
{
	double *vec = new double[length];
	return vec;
}

//------------------------------------------------------------------------
void vtkExtMisc::DeleteVector(double *v) 
{
	delete [] v;
}

//-------------------------------------------------------------------------
void vtkExtMisc::WriteMatrixASCII(double **m, int rows, int cols, const std::string& fname)
{
	std::ofstream out(fname.c_str(), std::ios::out | std::ios::trunc);

	if (out) 
	{
//		out << rows << " " << cols << std::endl;

		for (int r = 0; r < rows; r++) 
		{
			for (int c = 0; c < cols; c++) 
			{
				out << m[r][c] << " ";
			}
			out << std::endl;
		}
	} 
	else 
	{
		throw ("WriteMatrixASCII: could not open matrix file for writing");
	}
}


void vtkExtMisc::WriteMatrixVTKFormat(double **m, int rows, int cols, const std::string& fname, bool binary)
{

	vtkDenseArray<double>* matrix = vtkDenseArray<double>::New();  
	matrix->Resize(rows, cols);

	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			double val = m[r][c];
			matrix->SetValue(r, c, val);
		}
	}

	vtkArrayWriter::Write(matrix, fname, binary);

	matrix->Delete();
}

double ** vtkExtMisc::ReadMatrixVTKFormat(int &rows, int &cols, const std::string& fname)
{

	std::ifstream fist;

	// Has to been opened in binary
	fist.open(fname.c_str(), std::ifstream::binary);
//	std::ifstream fist(fname.c_str());
	if (!fist)
	{
		std::cerr << "Could not read " << std::endl;
		return NULL;
	}

	vtkArray *a = vtkArrayReader::Read(fist);
	if (!a)
	{
		std::cerr << "Could not read " << fname << std::endl;
		return NULL;
	}
	vtkDenseArray<double>* matrix = vtkDenseArray<double>::SafeDownCast(a);
	rows = matrix->GetExtent(0).GetSize();
	cols = matrix->GetExtent(1).GetSize();
	
	double **m = NewMatrix(rows, cols);
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			double val = matrix->GetValue(r, c);
			m[r][c] = val;
		}
	}
	a->Delete();
	return m;
}


void vtkExtMisc::WriteFloatArray(vtkFloatArray *a, const std::string &fname)
{
	std::ofstream fost(fname.c_str());
	for (int i = 0; i < a->GetNumberOfTuples(); i++)
	{
		fost << a->GetValue(i) << std::endl;
	}
}

bool vtkExtMisc::ReadDoubleArray(const std::string& fname, vtkDoubleArray *array)
{
	std::ifstream fist(fname.c_str());
	if (!fist)
	{
		std::cerr << "Could not read " << fname << std::endl;
		return false;
	}

	bool stop = false;
	do
	{
		std::string t1;
		std::getline(fist, t1);

		if (!fist.fail() && !fist.eof() && t1.size() > 0)
		{
			double val = atof(t1.c_str());
			array->InsertNextValue(val);
		}
		else
			stop = true;
	} while (!stop);
	if (array->GetNumberOfTuples() < 1)
	{
		std::cerr << "Less than 1 value read" << std::endl;
		return false;
	}
	return true;
}

void vtkExtMisc::TransformMatrixFromLine(double *p1, double *p2, double theta, vtkMatrix4x4* mat)
{
	double v1[3]; v1[0] = p2[0] - p1[0]; v1[1] = p2[1] - p1[1]; v1[2] = p2[2] - p1[2]; 

	vtkMath::Normalize(v1);
	
	double v2[3];
	double v3[3];

	vtkMath::Perpendiculars(v1, v2, v3, theta);

	mat->Identity();

	for (int i = 0; i < 3; i++)
	{
		mat->SetElement(i, 0, v1[i]);
		mat->SetElement(i, 1, v2[i]);
		mat->SetElement(i, 2, v3[i]);

		// Translation
		mat->SetElement(i, 3, p1[i]);
	}
}

void vtkExtMisc::WritePointsAsPCAFittedBlob(vtkPoints *pts, const std::string &fname)
{
	vtkPrincipalAxisTransform *princ = vtkPrincipalAxisTransform::New();
	princ->SetSource(pts);

	vtkSphereSource *sphere = vtkSphereSource::New();
	sphere->SetRadius(1);
	sphere->SetPhiResolution(20);
	sphere->SetThetaResolution(20);

	vtkTransformPolyDataFilter *tf = vtkTransformPolyDataFilter::New();
	tf->SetTransform(princ);
	tf->SetInputConnection(sphere->GetOutputPort());
	tf->Update();

	WritePDVTK(tf->GetOutput(), fname);

	tf->Delete();
	sphere->Delete();
	princ->Delete();
}

void vtkExtMisc::FitPlaneToPointsLSQ(vtkPoints *pts, vtkPlaneSource *plane, double size)
{
	// Taken from vtkPrincipalAxisTransform
	bool DoRotate = true;
	bool DoTranslate = true;
	bool DoScale = true;
	int i,j;

	vtkMatrix4x4 *Matrix = vtkMatrix4x4::New();
	
	if (pts->GetNumberOfPoints() == 0)
	{
		std:: cerr << "Can't execute with NULL or empty input";
		return;
	}
	
	//
	// Compute mean
	//
	vtkIdType numPts = pts->GetNumberOfPoints();

	vtkIdType pointId;

	double mean[3] = {0.0, 0.0, 0.0};
	for (pointId = 0; pointId < numPts; pointId++)
	{
		double x[3];
		pts->GetPoint(pointId, x);
		for (i = 0; i < 3; i++)
		{
			mean[i] += x[i];
		}
	}
	for (i=0; i < 3; i++)
	{
		mean[i] /= numPts;
	}
	

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
	for (pointId = 0; pointId < numPts; pointId++)
	{
		double x[3];
		pts->GetPoint(pointId, x);

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
	double *v[3], v0[3], v1[3], v2[3];
	double evals[3] = {0.0, 0.0, 0.0};
	// v is a vector of pointers pointing into the three vectors v0, v1 and v2
	// so basically v is a matrix with rows v0, v1 and v2
	// the eigenvectors are returned in v. the eigenvalues in evals;

	v[0] = v0; v[1] = v1; v[2] = v2; 
	vtkMath::Jacobi(a,evals,v);

//	cout << "Evals: " << evals[0] << " " << evals[1] << " " << evals[2] << endl;

	// Extract eigenvectors
	
	Matrix->Identity();
	if (DoRotate)
	{
		if (DoScale)
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					Matrix->SetElement(i, j, v[i][j] * sqrt(evals[j]));
				}
			}
		}
		else
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					Matrix->SetElement(i, j, v[i][j]);
				}
			}
		}
	}
	else if (DoScale)
	{
		for (i = 0; i < 3; i++)
		{
			Matrix->SetElement(i, i, sqrt(evals[i]));
		}
	}

	if (DoTranslate)
	{
		// Hard code translation
		Matrix->SetElement(0, 3, mean[0]);
		Matrix->SetElement(1, 3, mean[1]);
		Matrix->SetElement(2, 3, mean[2]);
	}

	// Check
	double (*matrix)[4] = Matrix->Element;
	double U[3][3], VT[3][3];
	
	for (i = 0; i < 3; i++) 
	{
		U[0][i] = matrix[0][i];
		U[1][i] = matrix[1][i];
		U[2][i] = matrix[2][i];
	}
	
	double scale[3];
	vtkMath::SingularValueDecomposition3x3(U, U, scale, VT);

	// This is dirty, but I do not know how else to avoid the scale being negative!!!
	if (scale[0] < 0 || scale[1] < 0 || scale[2] < 0)
	{
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				Matrix->SetElement(i, j, -(Matrix->GetElement(i,j)));
			}
		}
	}


	// Always set origin first
	plane->SetOrigin(mean[0], mean[1], mean[2]);
	plane->SetPoint1(mean[0] + Matrix->GetElement(0,1)*size, 
					 mean[1] + Matrix->GetElement(1,1)*size, 
					 mean[2] + Matrix->GetElement(2,1)*size);
	plane->SetPoint2(mean[0] + Matrix->GetElement(0,0)*size, 
					 mean[1] + Matrix->GetElement(1,0)*size, 
					 mean[2] + Matrix->GetElement(2,0)*size);

	// Only used for moving the plane
	plane->SetCenter(mean[0], mean[1], mean[2]);

	Matrix->Delete();
}
void vtkExtMisc::IterativeCovarianceMatrix(vtkPoints *pts, 
							   double *prevMean, double **prevCov, int prevNumberOfPoints,
							   double *newMean, double **newCov, int *numberOfPoints)
{
	int i;

	// With zeros in prev data the algo reduces to ordinary covariance calculation
	// but a bit too much is calculated in this way
	if(prevMean == NULL || prevCov == NULL)
	{
		std:: cerr << "previous mean or covariance is NULL pointer ";
		return;
	}

	if (pts->GetNumberOfPoints() == 0)
	{
		std:: cerr << "Can't execute with NULL or empty input";
		return;
	}
	//
	// Compute mean for new points
	//
	vtkIdType numPts = pts->GetNumberOfPoints();
	double mean[3] = {0.0, 0.0, 0.0};
	double x[3];
	for (vtkIdType pointId = 0; pointId < numPts; pointId++)
	{
		pts->GetPoint(pointId, x);
		for (i = 0; i < 3; i++)
		{
			mean[i] += x[i];
		}
	}
	for (i=0; i < 3; i++)
	{
		mean[i] /= numPts;
	}

	// Total number of points
	*numberOfPoints = numPts+prevNumberOfPoints;

	// ratios of new and previous points.
	double wnew = (double) numPts/(double)(*numberOfPoints);
	double wprev = (double) prevNumberOfPoints/(double)(*numberOfPoints);
	// new mean over all points
	for (i = 0; i < 3; i++)
	{
		newMean[i] = wnew*mean[i] + wprev *prevMean[i];
	}

	std::cout << "weight new and prev:" << wnew << "," << wprev << std::endl;
	std::cout << "old Mean in sub :" << prevMean[0] << "," << prevMean[1]<< "," << prevMean[2] << std::endl;
	std::cout << "Mean for new data in sub :" << mean[0] << "," << mean[1]<< "," << mean[2] << std::endl;
	std::cout << "Mean for all in sub :" << newMean[0] << "," << newMean[1]<< "," << newMean[2] << std::endl;

	// The mixed mean term weighted covariance of means
	double tmp = prevNumberOfPoints*wnew*wnew;
	double prevMinusCurMean[3];
	for (i=0; i < 3; i++)
	{
		prevMinusCurMean[i]=prevMean[i]-mean[i];
	}
	// MixedTerm 
	double mt[3][3];
	for(i=0; i<3; i++)
	{
		for(int j=0;j < 3; j++)
		{
			mt[i][j]=tmp*prevMinusCurMean[i]*prevMinusCurMean[j];
		}
	}


	// a is a vector of pointers pointing into the three vectors a0, a1 and a2
	// so basically a is a matrix with rows a0, a1 and a2
	double a[3][3] = {0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0, 0.0, 0.0};

	// Compute covariance matrix for new pts 
	// but don't use the mean over all points
	double xp[3];
	for (vtkIdType pointId = 0; pointId < numPts; pointId++)
	{
		pts->GetPoint(pointId, x);
		xp[0] = x[0] - newMean[0]; xp[1] = x[1] - newMean[1]; xp[2] = x[2] - newMean[2];
		for (i=0; i < 3; i++)
		{
			a[0][i] += xp[0] * xp[i];
			a[1][i] += xp[1] * xp[i];
			a[2][i] += xp[2] * xp[i];
		}
	}//for all points

  	for (i=0; i < 3; i++)
	{
		for(int j=0; j < 3; j++)
		{
			newCov[i][j] = a[i][j]+mt[i][j]+prevNumberOfPoints*prevCov[i][j];
		}
	}

	for (i=0; i < 3; i++)
	{
		for(int j=0; j < 3; j++)
		{
			newCov[i][j] /= (*numberOfPoints);
		}
	}

}
void vtkExtMisc::ComputeEigenVectorsOfPointCloud(vtkPoints *pts, double* mean, 
												 double *v0, double *v1, double *v2, double* evals)
{
	int i;

	vtkMatrix4x4 *Matrix = vtkMatrix4x4::New();
	
	if (pts->GetNumberOfPoints() == 0)
	{
		std:: cerr << "Can't execute with NULL or empty input";
		return;
	}
	
	//
	// Compute mean
	//
	vtkIdType numPts = pts->GetNumberOfPoints();

//	double mean[3] = {0.0, 0.0, 0.0};
	mean[0] = 0; mean[1] = 0; mean[2] = 0;
	for (vtkIdType pointId = 0; pointId < numPts; pointId++)
	{
		double x[3];
		pts->GetPoint(pointId, x);
		for (i = 0; i < 3; i++)
		{
			mean[i] += x[i];
		}
	}
	for (i=0; i < 3; i++)
	{
		mean[i] /= numPts;
	}
	

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
	for (vtkIdType pointId = 0; pointId < numPts; pointId++)
	{
		double x[3];
		pts->GetPoint(pointId, x);

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

// Adapted from https://math.stackexchange.com/questions/61719/finding-the-intersection-point-of-many-lines-in-3d-point-closest-to-all-lines
void vtkExtMisc::ComputeLeastSquaresIntersectionPoint(vtkPoints *startPoints, vtkPoints *endPoints, double *p)
{
	int nPoints = startPoints->GetNumberOfPoints();
	// Line normals
	std::vector<double> nx;
	std::vector<double> ny;
	std::vector<double> nz;
	for (int i = 0; i < nPoints; i++)
	{
		double dx = endPoints->GetPoint(i)[0] - startPoints->GetPoint(i)[0];
		double dy = endPoints->GetPoint(i)[1] - startPoints->GetPoint(i)[1];
		double dz = endPoints->GetPoint(i)[2] - startPoints->GetPoint(i)[2];

		double l = std::sqrt(dx * dx + dy * dy + dz * dz);

		if (l > 0)
		{
			nx.push_back(dx / l);
			ny.push_back(dy / l);
			nz.push_back(dz / l);
		}
	}
	size_t nValid = nx.size();
	if (nValid != nPoints)
	{
		// TODO: handle non valid lines
		std::cerr << "Line with no length given" << std::endl;
		return;
	}

	//SXX = sum(nx. ^ 2 - 1);
	//SYY = sum(ny. ^ 2 - 1);
	//SZZ = sum(nz. ^ 2 - 1);
	//SXY = sum(nx.*ny);
	//SXZ = sum(nx.*nz);
	//SYZ = sum(ny.*nz);

	double SXX = 0;
	double SYY = 0; 
	double SZZ = 0; 
	double SXY = 0;
	double SXZ = 0;
	double SYZ = 0;
	double CX = 0;
	double CY = 0;
	double CZ = 0;
	for (int i = 0; i < nValid; i++)
	{
		SXX += nx[i] * nx[i] - 1;
		SYY += ny[i] * ny[i] - 1;
		SZZ += nz[i] * nz[i] - 1;
		SXY += nx[i] * ny[i];
		SXZ += nx[i] * nz[i];
		SYZ += ny[i] * nz[i];
	
		double pa[3];
		startPoints->GetPoint(i, pa);

		CX += pa[0] * (nx[i] * nx[i] - 1) + pa[1] * nx[i] * ny[i]       + pa[2] * nx[i] * nz[i];
		CY += pa[0] * nx[i] * ny[i]       + pa[1] * (ny[i] * ny[i] - 1) + pa[2] * ny[i] * nz[i];
		CZ += pa[0] * nx[i] * nz[i]       + pa[1] * ny[i] * nz[i]       + pa[2] * (nz[i] * nz[i] - 1);
	}

	double **S = NewMatrix(3, 3);
	S[0][0] = SXX;  S[0][1] = SXY;  S[0][2] = SXZ;
	S[1][0] = SXY;  S[1][1] = SYY;  S[1][2] = SYZ;
	S[2][0] = SXZ;  S[2][1] = SYZ;  S[2][2] = SZZ;

	double **C = NewMatrix(3, 1);
	C[0][0] = CX;
	C[1][0] = CY;
	C[2][0] = CZ;

	double **PS = NewMatrix(3, 1);

	vtkMath::SolveLeastSquares(3, S, 3, C, 1, PS);
	p[0] = PS[0][0];
	p[1] = PS[1][0];
	p[2] = PS[2][0];
}

/* 
while iterations < k {
maybeinliers = n randomly selected values from data
maybemodel = model parameters fitted to maybeinliers
alsoinliers = empty set
for every point in data not in maybeinliers {
if point fits maybemodel with an error smaller than t
add point to alsoinliers
}
if the number of elements in alsoinliers is > d {
% this implies that we may have found a good model
% now test how good it is
bettermodel = model parameters fitted to all points in maybeinliers and alsoinliers
thiserr = a measure of how well model fits these points
if thiserr < besterr {
bestfit = bettermodel
besterr = thiserr
}
}
increment iterations
}
*/

void vtkExtMisc::ComputeLeastSquaresIntersectionPointWithRANSAC(vtkPoints *startPoints, vtkPoints *endPoints, double *p)
{
	int iterations = 100;
	double bestError = 10000000000;
	double bestP[3];
	double distThres = 20;
	int d = 10;

	int nLines = startPoints->GetNumberOfPoints();
	d = nLines / 2;
	bestP[0] = 0;
	bestP[1] = 0; 
	bestP[2] = 0;
	int usedLines = 0;

	for (int i = 0; i < iterations; i++)
	{
		// Maybeinliers
		int n1 = (int)(vtkMath::Random() * nLines);
		int n2 = 0;
		do 
		{
			n2 = (int)(vtkMath::Random() * nLines);
		} while (n1 == n2);
		vtkPoints *SPt = vtkPoints::New();
		vtkPoints *EPt = vtkPoints::New();

		SPt->InsertNextPoint(startPoints->GetPoint(n1));
		SPt->InsertNextPoint(startPoints->GetPoint(n2));
		EPt->InsertNextPoint(endPoints->GetPoint(n1));
		EPt->InsertNextPoint(endPoints->GetPoint(n2));

		// maybeinliers
		double Pest[3];
		ComputeLeastSquaresIntersectionPoint(SPt, EPt, Pest);

		std::vector<int> alsoinliers;
		for (int j = 0; j < nLines; j++)
		{
			if (j != n1  && j != n2)
			{
				double dist2 = vtkLine::DistanceToLine(Pest, startPoints->GetPoint(j), endPoints->GetPoint(j));
				if (dist2 < distThres * distThres)
				{
					alsoinliers.push_back(j);

					SPt->InsertNextPoint(startPoints->GetPoint(j));
					EPt->InsertNextPoint(endPoints->GetPoint(j));
				}
			}
		}
		if (alsoinliers.size() > d)
		{
			// Reestimae parameters
			double Pest[3];
			ComputeLeastSquaresIntersectionPoint(SPt, EPt, Pest);

			double sumSQ = 0;
			for (int j = 0; j < SPt->GetNumberOfPoints(); j++)
			{
				double dist2 = vtkLine::DistanceToLine(Pest, SPt->GetPoint(j), EPt->GetPoint(j));
				sumSQ += dist2;
			}
			sumSQ /= SPt->GetNumberOfPoints();
			if (sumSQ < bestError)
			{
				bestError = sumSQ;
				bestP[0] = Pest[0];
				bestP[1] = Pest[1];
				bestP[2] = Pest[2];
				usedLines = SPt->GetNumberOfPoints();
			}

		}
		SPt->Delete();
		EPt->Delete();
	}
	std::cout << "RANSAC best error " << bestError << " using " << usedLines << " lines out of " << startPoints->GetNumberOfPoints() << std::endl;
	p[0] = bestP[0];
	p[1] = bestP[1];
	p[2] = bestP[2];
}

void vtkExtMisc::SavePlaneAsCorners(vtkPlaneSource *plane, const std::string& fname)
{
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkPolyData *pd = vtkPolyData::New();

	vtkIdType id = pts->InsertNextPoint(plane->GetOrigin());
	verts->InsertNextCell(1);
	verts->InsertCellPoint(id);
	
	id = pts->InsertNextPoint(plane->GetPoint1());
	verts->InsertNextCell(1);
	verts->InsertCellPoint(id);

	id = pts->InsertNextPoint(plane->GetPoint2());
	verts->InsertNextCell(1);
	verts->InsertCellPoint(id);

	pd->SetPoints(pts);
	pd->SetVerts(verts);

	vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
	pWriter->SetInputData(pd);
	pWriter->SetFileName(fname.c_str());
	pWriter->Write();

	pWriter->Delete();
	pts->Delete();
	verts->Delete();
	pd->Delete();
}

void vtkExtMisc::Plane2Parameters(double *po, double *p1, double *p2, double &SideLength, double &A, double &B, double &C, double &D)
{
	vtkPlaneSource *plane = vtkPlaneSource::New();
	 plane->SetOrigin(po);
	 plane->SetPoint1(p1);
	 plane->SetPoint2(p2);
	 plane->Update();

	SideLength = (sqrt(vtkMath::Distance2BetweenPoints(po, p1)) + sqrt(vtkMath::Distance2BetweenPoints(po, p2))) / 2;

	double  pNorm[3];
	plane->GetNormal(pNorm);
	A = pNorm[0];
	B = pNorm[1];
	C = pNorm[2];
	D = -(A * po[0] + B * po[1] + C * po[2]);

	// Normalise it
	double Plength = sqrt(A*A + B*B + C*C + D*D);
	A /= Plength;B /= Plength;C /= Plength;	D /= Plength;

	plane->Delete();
}

void vtkExtMisc::FindVertexNeighbours(vtkPolyData* mesh, int VertexID, vtkIdList *NeighbourIDs, double NeighborhoodSize2)
{
	double PC[3];
	mesh->GetPoint(VertexID, PC);

    vtkIdList *verts  = vtkIdList::New();
	vtkIdList *cellids = vtkIdList::New();

	std::deque<vtkIdType> visitlist;

	visitlist.push_back(VertexID);

	// 0 = not visited
	// 1 = in que
	// 2 = processed
	std::vector<int> statusVec(mesh->GetNumberOfPoints(), 0);

	while (!visitlist.empty())
	{
		vtkIdType curId = visitlist.front(); 
		visitlist.pop_front();
		
		statusVec[curId] = 2;

		// Get the cells that point with curID belongs to:
		mesh->GetPointCells(curId, cellids);
		for (int nc = 0; nc < cellids->GetNumberOfIds(); nc++)
		{
			// Cell ID
			vtkIdType cid = cellids->GetId(nc);

			// Get the points of the neighbour cell
			// it is assumed that all points in a cell is connected...
			mesh->GetCellPoints(cid, verts);

			// Iterate over all neigbhour points
			for (int nv = 0; nv < verts->GetNumberOfIds(); nv++)
			{
				vtkIdType tid = verts->GetId(nv);

				if (statusVec[tid] == 0)
				{
					statusVec[tid] = 1;

					double Pcur[3];
					mesh->GetPoint(tid, Pcur);

					if (vtkMath::Distance2BetweenPoints(PC, Pcur) < NeighborhoodSize2)
					{
						visitlist.push_back(tid);
						NeighbourIDs->InsertNextId(tid);
					}
				}
			}
		}
	}
	verts->Delete();
	cellids->Delete();
}


vtkPolyDataReader *vtkExtMisc::SafeReadPolyData(const std::string& fname)
{
	vtkPolyDataReader *read = vtkPolyDataReader::New();
	 read->SetFileName(fname.c_str());
	 read->Update();

	if (read->GetOutput()->GetNumberOfPoints() == 0)
	{
		std::cerr << "Could not read: " << fname << std::endl;
		read->Delete();
		read = NULL;
	}
	return read;
}

vtkSTLReader *vtkExtMisc::SafeReadPolyDataSTL(const std::string& fname)
{
	vtkSTLReader *read = vtkSTLReader::New();
	 read->SetFileName(fname.c_str());
	 read->Update();

	if (read->GetOutput()->GetNumberOfPoints() == 0)
	{
		std::cerr << "Could not read: " << fname << std::endl;
		read->Delete();
		read = NULL;
	}
	return read;
}

vtkOBJReader *vtkExtMisc::SafeReadPolyDataOBJ(const std::string& fname)
{
	vtkOBJReader *read = vtkOBJReader::New();
	 read->SetFileName(fname.c_str());
	 read->Update();

	if (read->GetOutput()->GetNumberOfPoints() == 0)
	{
		std::cerr << "Could not read: " << fname << std::endl;
		read->Delete();
		read = NULL;
	}
	return read;
}


void vtkExtMisc::CenterOfMass(vtkPolyData* pd, vtkIdList* pids, double *CM)
{
	CM[0] = 0; CM[1] = 0; CM[2] = 0;

	const int N = pids->GetNumberOfIds();

	if (!N) return;

	for (int i = 0; i < N; i++)
	{
		double pNi[3];
		pd->GetPoint(pids->GetId(i), pNi);

		CM[0] += pNi[0];
		CM[1] += pNi[1];
		CM[2] += pNi[2];
	}
	CM[0] /= N; CM[1] /= N; CM[2] /= N;
}

void vtkExtMisc::CenterOfMass(vtkPolyData* pd, double *CM)
{
	CM[0] = 0; CM[1] = 0; CM[2] = 0;

	const int N = pd->GetNumberOfPoints();

	if (!N) return;

	for (int i = 0; i < N; i++)
	{
		double pNi[3];
		pd->GetPoint(i, pNi);

		CM[0] += pNi[0];
		CM[1] += pNi[1];
		CM[2] += pNi[2];
	}
	CM[0] /= N; CM[1] /= N; CM[2] /= N;
}

void vtkExtMisc::CenterOfMass(vtkPoints* pd, double *CM)
{
	CM[0] = 0; CM[1] = 0; CM[2] = 0;

	const int N = pd->GetNumberOfPoints();

	if (!N) return;

	for (int i = 0; i < N; i++)
	{
		double pNi[3];
		pd->GetPoint(i, pNi);

		CM[0] += pNi[0];
		CM[1] += pNi[1];
		CM[2] += pNi[2];
	}
	CM[0] /= N; CM[1] /= N; CM[2] /= N;
}

void vtkExtMisc::TranslateCenterOfMassToOrigo(vtkPolyData *inPD, vtkPolyData *outPD)
{
	double CM[3];
	vtkExtMisc::CenterOfMass(inPD, CM);

	vtkTransform *trans = vtkTransform::New();
		trans->Translate(-CM[0], -CM[1], -CM[2]);

	vtkTransformPolyDataFilter *tfilt = vtkTransformPolyDataFilter::New();
	tfilt->SetInputData(inPD);
	tfilt->SetTransform(trans);
	tfilt->Update();

	outPD->DeepCopy(tfilt->GetOutput());

	trans->Delete();
	tfilt->Delete();
}


void vtkExtMisc::NormaliseSize( vtkPolyData *inPD, vtkPolyData *outPD, double size, double *CM, double &scale )
{
	vtkExtMisc::CenterOfMass(inPD, CM);

	double bounds[6];
	inPD->GetBounds(bounds);

	double maxSide = std::max(bounds[1]-bounds[0], std::max(bounds[3]-bounds[2], bounds[5]-bounds[4]));
	scale = size / maxSide;

	vtkTransform *trans = vtkTransform::New();
	trans->Scale(scale, scale,scale);
	trans->Translate(-CM[0], -CM[1], -CM[2]);

	vtkTransformPolyDataFilter *tfilt = vtkTransformPolyDataFilter::New();
	tfilt->SetInputData(inPD);
	tfilt->SetTransform(trans);
	tfilt->Update();

	outPD->DeepCopy(tfilt->GetOutput());

	trans->Delete();
	tfilt->Delete();

}

/*
http://www.exaflop.org/docs/cgafaq/cga2.html:

To find the area of a planar polygon not in the x-y plane, use:

2 A(P) = abs(N . (sum_{i=0}^{n-1} (v_i x v_{i+1})))
where N is a unit vector normal to the plane. The `.' represents the dot product operator, the `x' represents the cross product operator, and abs() is the absolute value function. 

[Gems II] pp. 170-171:

"Area of Planar Polygons and Volume of Polyhedra", Ronald N. Goldman
*/

double vtkExtMisc::PolygonArea(vtkPoints *poly)
{
	const int N = poly->GetNumberOfPoints();

	if (N < 3)
		return 0;

	double sumVec[] = {0,0,0};

	// Make the sum of crossproducts
	for (int i = 0; i < N-1; i++)
	{
		double tempVec[3];

		double pNi[3];
		poly->GetPoint(i, pNi);
		vtkMath::Cross(pNi, poly->GetPoint(i+1), tempVec);

		sumVec[0] += tempVec[0]; sumVec[1] += tempVec[1]; sumVec[2] += tempVec[2]; 
	}

	// Calculate a plane normal
	double v1[3];
	double v2[3];

	int N1 = N/3;
	int N2 = 2*N/3;
	if (N1 == N2)
	{
		N1 = 1;
		N2 = 2;
	}


	double pN1[3];
	poly->GetPoint(N1, pN1);
	double p0[3];
	poly->GetPoint(0, p0);
	double pN2[3];
	poly->GetPoint(N2, pN2);

	// we guess that the vector from point 0 to point 1 is not
	// parallel to the vector from point 0 to point N/2
	v1[0] = pN1[0] - p0[0];
	v1[1] = pN1[1] - p0[1];
	v1[2] = pN1[2] - p0[2];
	
	v2[0] = pN2[0] - p0[0];
	v2[1] = pN2[1] - p0[1];
	v2[2] = pN2[2] - p0[2];

	double PolyNormal[3];
	vtkMath::Cross(v1, v2, PolyNormal);
	vtkMath::Normalize(PolyNormal);
	
	double Area = fabs(vtkMath::Dot(sumVec, PolyNormal)) / 2;

	return Area;
}


// http://geomalgorithms.com/a01-_area.html
double vtkExtMisc::PolygonAreaByProjection(vtkPoints *poly)
{
	int i,j,k;

	const int N = poly->GetNumberOfPoints();

	if (N < 3)
		return 0;

	// Calculate a plane normal
	double v1[3];
	double v2[3];

	int N1 = N/3;
	int N2 = 2*N/3;
	if (N1 == N2)
	{
		N1 = 1;
		N2 = 2;
	}

	double pN1[3];
	poly->GetPoint(N1, pN1);
	double p0[3];
	poly->GetPoint(0, p0);
	double pN2[3];
	poly->GetPoint(N2, pN2);

	// we guess that the vector from point 0 to point 1 is not
	// parallel to the vector from point 0 to point N/2
	v1[0] = pN1[0] - p0[0];
	v1[1] = pN1[1] - p0[1];
	v1[2] = pN1[2] - p0[2];
	
	v2[0] = pN2[0] - p0[0];
	v2[1] = pN2[1] - p0[1];
	v2[2] = pN2[2] - p0[2];

	double PolyNormal[3];
	vtkMath::Cross(v1, v2, PolyNormal);
	vtkMath::Normalize(PolyNormal);

    // select largest abs coordinate to ignore for projection
    double ax = (PolyNormal[0] > 0 ? PolyNormal[0] : -PolyNormal[0]);     // abs x-coord
    double ay = (PolyNormal[1] > 0 ? PolyNormal[1] : -PolyNormal[1]);     // abs y-coord
    double az = (PolyNormal[2] > 0 ? PolyNormal[2] : -PolyNormal[2]);     // abs z-coord

	// coord to ignore: 1=x, 2=y, 3=z
    int coord = 3;
    if (ax > ay) 
	{
        if (ax > az) 
			coord = 1;    // ignore x-coord
    }
    else if (ay > az) 
		coord = 2;   // ignore y-coord


	double Area = 0;

    // compute area of the 2D projection
    for (i=1, j=2, k=0; i <= N-2; i++, j++, k++)
	{
		double pNi[3];
		double pNj[3];
		double pNk[3];
		poly->GetPoint(i, pNi);
		poly->GetPoint(j, pNj);
		poly->GetPoint(k, pNk);
        switch (coord) 
		{
			case 1:
				Area += (pNi[1] * (pNj[2] - pNk[2]));
				continue;
			case 2:
				Area += (pNi[0] * (pNj[2] - pNk[2]));
				continue;
			case 3:
				Area += (pNi[0] * (pNj[1] - pNk[1]));
				continue;
		}
    }

    // scale to get area before projection
    double an = sqrt( ax*ax + ay*ay + az*az);  // length of normal vector
    switch (coord) 
	{
    case 1:
        Area *= (an / (2*ax));
        break;
    case 2:
        Area *= (an / (2*ay));
        break;
    case 3:
        Area *= (an / (2*az));
    }

	return fabs(Area);
}
/*
// http://geometryalgorithms.com/Archive/algorithm_0101/#3D%20Polygons
double vtkExtMisc::PolygonAreaByProjection(vtkCell *poly)
{
	int i,j,k;

	const int N = poly->GetNumberOfPoints();

	// Look-up-table to points
	vtkIdList *pLUT = poly->GetPointIds();

	if (N < 3)
		return 0;

	// Calculate a plane normal
	double v1[3];
	double v2[3];

	int N1 = N/3;
	int N2 = 2*N/3;
	if (N1 == N2)
	{
		N1 = 1;
		N2 = 2;
	}

	vtkPoints *pts = poly->GetPoints();

	pts->GetPoint(pLUT->GetId(N1));

	// we guess that the vector from point 0 to point 1 is not
	// parallel to the vector from point 0 to point N/2
	v1[0] = pts->GetPoint(pLUT->GetId(N1))[0] - pts->GetPoint(pLUT->GetId(0))[0];
	v1[1] = pts->GetPoint(pLUT->GetId(N1))[1] - pts->GetPoint(pLUT->GetId(0))[1];
	v1[2] = pts->GetPoint(pLUT->GetId(N1))[2] - pts->GetPoint(pLUT->GetId(0))[2];
	
	v2[0] = pts->GetPoint(pLUT->GetId(N2))[0] - pts->GetPoint(pLUT->GetId(0))[0];
	v2[1] = pts->GetPoint(pLUT->GetId(N2))[1] - pts->GetPoint(pLUT->GetId(0))[1];
	v2[2] = pts->GetPoint(pLUT->GetId(N2))[2] - pts->GetPoint(pLUT->GetId(0))[2];

	double PolyNormal[3];
	vtkMath::Cross(v1, v2, PolyNormal);
	vtkMath::Normalize(PolyNormal);

    // select largest abs coordinate to ignore for projection
    double ax = (PolyNormal[0] > 0 ? PolyNormal[0] : -PolyNormal[0]);     // abs x-coord
    double ay = (PolyNormal[1] > 0 ? PolyNormal[1] : -PolyNormal[1]);     // abs y-coord
    double az = (PolyNormal[2] > 0 ? PolyNormal[2] : -PolyNormal[2]);     // abs z-coord

	// coord to ignore: 1=x, 2=y, 3=z
    int coord = 3;
    if (ax > ay) 
	{
        if (ax > az) 
			coord = 1;    // ignore x-coord
    }
    else if (ay > az) 
		coord = 2;   // ignore y-coord


	double Area = 0;

    // compute area of the 2D projection
    for (i=1, j=2, k=0; i <= N-2; i++, j++, k++)
	{
        switch (coord) 
		{
			case 1:
				Area += (pts->GetPoint(pLUT->GetId(i))[1] * (pts->GetPoint(pLUT->GetId(j))[2] - pts->GetPoint(pLUT->GetId(k))[2]));
				continue;
			case 2:
				Area += (pts->GetPoint(pLUT->GetId(i))[0] * (pts->GetPoint(pLUT->GetId(j))[2] - pts->GetPoint(pLUT->GetId(k))[2]));
				continue;
			case 3:
				Area += (pts->GetPoint(pLUT->GetId(i))[0] * (pts->GetPoint(pLUT->GetId(j))[1] - pts->GetPoint(pLUT->GetId(k))[1]));
				continue;
		}
    }

    // scale to get area before projection
    double an = sqrt( ax*ax + ay*ay + az*az);  // length of normal vector
    switch (coord) 
	{
    case 1:
        Area *= (an / (2*ax));
        break;
    case 2:
        Area *= (an / (2*ay));
        break;
    case 3:
        Area *= (an / (2*az));
    }

	return Area;
}
*/
/*
// http://geometryalgorithms.com/Archive/algorithm_0101/#3D%20Polygons
double vtkExtMisc::PolygonAreaByProjection(vtkPoints *poly, double *Normal)
{
	int i,j,k;

	const int N = poly->GetNumberOfPoints();

	if (N < 3)
		return 0;

	// Calculate a plane normal
	double v1[3];
	double v2[3];

	int N1 = N/3;
	int N2 = 2*N/3;
	if (N1 == N2)
	{
		N1 = 1;
		N2 = 2;
	}

	// we guess that the vector from point 0 to point 1 is not
	// parallel to the vector from point 0 to point N/2
	v1[0] = pN1[0] - p0[0];
	v1[1] = pN1[1] - p0[1];
	v1[2] = pN1[2] - p0[2];
	
	v2[0] = pN2[0] - p0[0];
	v2[1] = pN2[1] - p0[1];
	v2[2] = pN2[2] - f[2];

	double PolyNormal[3];
	vtkMath::Cross(v1, v2, PolyNormal);
	vtkMath::Normalize(PolyNormal);

	double CheckNormal = vtkMath::Dot(Normal, PolyNormal);
	std::cout << "Normals check: " << CheckNormal << std::endl;


    // select largest abs coordinate to ignore for projection
    double ax = (PolyNormal[0] > 0 ? PolyNormal[0] : -PolyNormal[0]);     // abs x-coord
    double ay = (PolyNormal[1] > 0 ? PolyNormal[1] : -PolyNormal[1]);     // abs y-coord
    double az = (PolyNormal[2] > 0 ? PolyNormal[2] : -PolyNormal[2]);     // abs z-coord

	// coord to ignore: 1=x, 2=y, 3=z
    int coord = 3;
    if (ax > ay) 
	{
        if (ax > az) 
			coord = 1;    // ignore x-coord
    }
    else if (ay > az) 
		coord = 2;   // ignore y-coord


	double Area = 0;

    // compute area of the 2D projection
    for (i=1, j=2, k=0; i <= N; i++, j++, k++)
	{
        switch (coord) 
		{
			case 1:
				Area += (pNi[1] * (pNj[2] - pNk[2]));
				continue;
			case 2:
				Area += (pNi[0] * (pNj[2] - pNk[2]));
				continue;
			case 3:
				Area += (pNi[0] * (pNj[1] - pNk[1]));
				continue;
		}
    }

    // scale to get area before projection
    double an = sqrt( ax*ax + ay*ay + az*az);  // length of normal vector
    switch (coord) 
	{
    case 1:
        Area *= (an / (2*ax));
        break;
    case 2:
        Area *= (an / (2*ay));
        break;
    case 3:
        Area *= (an / (2*az));
    }

	return Area;
}
*/
/*
// area3D_Polygon(): computes the area of a 3D planar polygon
//    Input:  int n = the number of vertices in the polygon
//            Point* V = an array of n+2 vertices in a plane
//                       with V[n]=V[0] and V[n+1]=V[1]
//            Point N = unit normal vector of the polygon's plane
//    Return: the (float) area of the polygon
float
area3D_Polygon( int n, Point* V, Point N )
{
    float area = 0;
    float an, ax, ay, az;  // abs value of normal and its coords
    int   coord;           // coord to ignore: 1=x, 2=y, 3=z
    int   i, j, k;         // loop indices

    // select largest abs coordinate to ignore for projection
    ax = (N.x>0 ? N.x : -N.x);     // abs x-coord
    ay = (N.y>0 ? N.y : -N.y);     // abs y-coord
    az = (N.z>0 ? N.z : -N.z);     // abs z-coord

    coord = 3;                     // ignore z-coord
    if (ax > ay) {
        if (ax > az) coord = 1;    // ignore x-coord
    }
    else if (ay > az) coord = 2;   // ignore y-coord

    // compute area of the 2D projection
    for (i=1, j=2, k=0; i<=n; i++, j++, k++)
        switch (coord) {
        case 1:
            area += (V[i].y * (V[j].z - V[k].z));
            continue;
        case 2:
            area += (V[i].x * (V[j].z - V[k].z));
            continue;
        case 3:
            area += (V[i].x * (V[j].y - V[k].y));
            continue;
        }

    // scale to get area before projection
    an = sqrt( ax*ax + ay*ay + az*az);  // length of normal vector
    switch (coord) {
    case 1:
        area *= (an / (2*ax));
        break;
    case 2:
        area *= (an / (2*ay));
        break;
    case 3:
        area *= (an / (2*az));
    }
    return area;
}

*/


void vtkExtMisc::DistanceFromCMStats(vtkPoints *poly, double &DistMean, double &DistSDev)
{
	int i;
	DistMean = 0;
	DistSDev = 0;

	double CM[3];
	CM[0] = 0; CM[1] = 0; CM[2] = 0;

	const int N = poly->GetNumberOfPoints();

	if (!N) return;

	for (i = 0; i < N; i++)
	{
		double pNi[3];
		poly->GetPoint(i, pNi);
		CM[0] += pNi[0];
		CM[1] += pNi[1];
		CM[2] += pNi[2];
	}
	CM[0] /= N; CM[1] /= N; CM[2] /= N;

	double S = 0;
	double SS = 0;
	for (i = 0; i < N; i++)
	{
		double pNi[3];
		poly->GetPoint(i, pNi);
		double D2 = vtkMath::Distance2BetweenPoints(pNi, CM);
		S += sqrt(D2);
		SS += D2;
	}
	DistMean = S/N;
	DistSDev = sqrt(SS/N - DistMean*DistMean);
}

void vtkExtMisc::ProjectPolygonTo2DWriteToFile(vtkPoints *poly, const std::string fname)
{
	int i;
	const int N = poly->GetNumberOfPoints();

	if (N < 3)
	{
		std::cerr << "Not enough points in poly. Did not write " << fname << std::endl;
		return;
	}

	// Calculate a plane normal
	double v1[3];
	double v2[3];

	int N1 = N/3;
	int N2 = 2*N/3;
	if (N1 == N2)
	{
		N1 = 1;
		N2 = 2;
	}

	double pN1[3];
	poly->GetPoint(N1, pN1);
	double p0[3];
	poly->GetPoint(0, p0);
	double pN2[3];
	poly->GetPoint(N2, pN2);

	// we guess that the vector from point 0 to point 1 is not
	// parallel to the vector from point 0 to point N/2
	v1[0] = pN1[0] - p0[0];
	v1[1] = pN1[1] - p0[1];
	v1[2] = pN1[2] - p0[2];
	
	v2[0] = pN2[0] - p0[0];
	v2[1] = pN2[1] - p0[1];
	v2[2] = pN2[2] - p0[2];

	vtkMath::Normalize(v1);
	vtkMath::Normalize(v2);

	double PolyNormal[3];
	vtkMath::Cross(v1, v2, PolyNormal);
	vtkMath::Normalize(PolyNormal);

	// Recalculate v2
	vtkMath::Cross(v1, PolyNormal,v2);
	vtkMath::Normalize(v2);

	double Org[3];
	Org[0] = 0; Org[1] = 0; Org[2] = 0;
	// Calculate origin
	for (i = 0; i < N; i++)
	{
		double pNi[3];
		poly->GetPoint(i, pNi);
		Org[0] += pNi[0]; 
		Org[1] += pNi[1];
		Org[2] += pNi[2];
	}
	Org[0] /= N; Org[1] /= N; Org[2] /= N;


	std::ofstream ost(fname.c_str());
	if (!ost)
	{
		std::cerr << "Could not write " << fname << std::endl;
		return;
	}

    // compute area of the 2D projection
    for (i=0; i < N; i++)
	{
		double P[3];
		poly->GetPoint(i, P);
		P[0] -= Org[0];
		P[1] -= Org[1];
		P[2] -= Org[2];
		double xn = vtkMath::Dot(v1, P);
		double yn = vtkMath::Dot(v2, P);
		ost << xn << ", " << yn << std::endl;
    }
}


bool vtkExtMisc::InterSectionBetweenPolyLineAndPlane(vtkPolyData* Pline, double *pNorm, double *pOrg,
																int& lowPointId, double &IntersectT, double *iPoint)
{
	int i;

	lowPointId = 0;
	IntersectT = 0;
	for (i = 0; i < Pline->GetNumberOfPoints()-1; i++)
	{
		double p1[3];
		Pline->GetPoint(i,p1);
		double p2[3];
		Pline->GetPoint(i+1,p2);

		double t = 0;
		int Intersect = vtkPlane::IntersectWithLine(p1, p2, pNorm, pOrg, t, iPoint);

		if (Intersect)
		{
			lowPointId = i;
			IntersectT = t;
			return true;
		}
	}
	return false;
}

bool vtkExtMisc::MoveDistanceAlongPolyLine(vtkPolyData* Pline, int lowPointId, double t, double dist,
											int& NewLowPointId, double &NewT, double *NewP)
{
	if (lowPointId+1 >= Pline->GetNumberOfPoints())
		return false;

	// Find original point
	double Plo[3];
	double Phi[3];
	Pline->GetPoint(lowPointId, Plo);
	Pline->GetPoint(lowPointId+1, Phi);

	double Por[3];
	Por[0] = Plo[0] + t * (Phi[0] - Plo[0]);
	Por[1] = Plo[1] + t * (Phi[1] - Plo[1]);
	Por[2] = Plo[2] + t * (Phi[2] - Plo[2]);


	// Case 1: New point is on the same line
	double d = sqrt(vtkMath::Distance2BetweenPoints(Por, Phi));

	if (d > dist)
	{
		double tinc = dist/sqrt(vtkMath::Distance2BetweenPoints(Plo,Phi));
		NewLowPointId = lowPointId;
		NewT = t+tinc;
		NewP[0] = Plo[0] + NewT * (Phi[0] - Plo[0]);
		NewP[1] = Plo[1] + NewT * (Phi[1] - Plo[1]);
		NewP[2] = Plo[2] + NewT * (Phi[2] - Plo[2]);
		return true;
	}
	dist -= d;

	// Now we are placed at the 1. "high point"
	lowPointId++;

	while (lowPointId+1 < Pline->GetNumberOfPoints())
	{
		Pline->GetPoint(lowPointId, Plo);
		Pline->GetPoint(lowPointId+1, Phi);

		// length of current segment
		d = sqrt(vtkMath::Distance2BetweenPoints(Plo, Phi));

		if (d > dist)
		{
			// point is in this segment
			NewT = dist/sqrt(vtkMath::Distance2BetweenPoints(Plo,Phi));
			NewLowPointId = lowPointId;
			NewP[0] = Plo[0] + NewT * (Phi[0] - Plo[0]);
			NewP[1] = Plo[1] + NewT * (Phi[1] - Plo[1]);
			NewP[2] = Plo[2] + NewT * (Phi[2] - Plo[2]);
			return true;
		}
		dist -= d;
		lowPointId++;
	}
	std::cerr << "Exceeded the end of the path. Distance remaining: " << dist << std::endl;
	return false;
}

double vtkExtMisc::PointInsideOrOutside(double *porg, double *pnormal, double *p)
{
	double v[3];
	v[0] = p[0] - porg[0];
	v[1] = p[1] - porg[1];
	v[2] = p[2] - porg[2];

	return vtkMath::Dot(v, pnormal);
}

double vtkExtMisc::CutAndCalculateArea(vtkPolyData *pd, vtkPlane *plane)
{
	int i;
	double area = 0;

	vtkCutter *cutter = vtkCutter::New();
	 cutter->SetInputData(pd);
	 cutter->SetCutFunction(plane);
	 cutter->GenerateCutScalarsOff();
	 
	vtkStripper *stripper = vtkStripper::New();
	 stripper->SetInputConnection(cutter->GetOutputPort());
	 stripper->Update();

	// Take care of cut that consists of several segments
	int MaxCutNumber = 0;
	for (i = 0; i < stripper->GetOutput()->GetNumberOfCells(); i++)
	{
		double tL = vtkExtMisc::PolygonAreaByProjection(stripper->GetOutput()->GetCell(i)->GetPoints());
		if (tL > area)
		{
			area = tL;
			MaxCutNumber = i;
		}
	}
	
	cutter->Delete();
	stripper->Delete();

	return area;
}

double vtkExtMisc::SurfaceArea(vtkPolyData *pd)
{
	double area = 0;

	int NP = pd->GetNumberOfCells();

	if (NP == 0)
	{
		return 0;
	}

	for (int i = 0; i < NP; i++)
	{
		vtkTriangle *tri = vtkTriangle::SafeDownCast(pd->GetCell(i));
		if (!tri)
		{
			std::cerr << "Cell " << i << " is not a triangle!" << std::endl;
		}
		else
		{
			
			int id0 = tri->GetPointId(0);
			int id1 = tri->GetPointId(1);
			int id2 = tri->GetPointId(2);
			double p0[3];
			double p1[3];
			double p2[3];
			pd->GetPoint(id0, p0);
			pd->GetPoint(id1, p1);
			pd->GetPoint(id2, p2);

			double a = vtkTriangle::TriangleArea(p0, p1, p2);;
			area += a;
		}
	}
	
	return area;
}

double vtkExtMisc::FindClosestPointOnSurface(vtkPolyData *pd, double *x, double *cp)
{

	int subId;
	double pcoords[3], point[3],  weights[6];
  
	const int N = pd->GetNumberOfCells();
	if (N < 1)
	{
		std::cerr << "No cells in polydata" << std::endl;
		return 0;
	}

	double mindist2 = 0;
    int stat = pd->GetCell(0)->EvaluatePosition(x, point, subId, pcoords,
                  mindist2, weights);

	for (int i = 0; i < N; i++)
	{
		double dist2 = 0;

		int stat = pd->GetCell(i)->EvaluatePosition(x, point, subId, pcoords,
					  dist2, weights);
		if (stat != -1 && dist2 < mindist2)
		{
			mindist2 = dist2;
		}
	}

	return sqrt(mindist2);
}


std::string vtkExtMisc::GetSurfaceValues(vtkPolyData *pd, const std::string& sep, const std::string& spa)
{
	std::ostringstream ost;
	if (pd->GetNumberOfPoints())
		ost << "Number of points: " << spa << pd->GetNumberOfPoints() << sep;
	else
		ost << "Strange surface -> No points!" << sep;
	if (pd->GetNumberOfPolys())
		ost << "Number of polys: " << spa << pd->GetNumberOfPolys() << sep;
	if (pd->GetNumberOfLines())
		ost << "Number of lines: " << spa << pd->GetNumberOfLines() << sep;
	if (pd->GetNumberOfVerts())
		ost << "Number of verts: " << spa << pd->GetNumberOfVerts() << sep;
	if (pd->GetNumberOfStrips())
		ost << "Number of strips: " << spa << pd->GetNumberOfStrips() << sep;
	if (pd->GetNumberOfPieces())
		ost << "Number of pieces: " << spa << pd->GetNumberOfPieces() << sep;

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetScalars());
	if (!scalars)
	{
		if (pd->GetPointData()->GetScalars() && pd->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
		{
			ost << "Point scalars: 3 component scalars - probably RGB" << sep;
		}
		else if (pd->GetPointData()->GetScalars())
		{
			ost << "Point scalars: " << pd->GetPointData()->GetScalars()->GetNumberOfComponents() << " component scalars" << sep;
		}
		else
			ost << "No point scalars" << sep;
	}
	else
	{
		std::vector<double> scals(pd->GetNumberOfPoints());

		for (int i = 0; i < pd->GetNumberOfPoints(); i++)
		{
			scals[i] = scalars->GetValue(i);
		}
		
		double frac05 = 0;
		double frac95 = 0;
		double median = 0;
		double mean = 0;
		double sdev = 0;
		double minScal = 0;
		double maxScal = 0;
		CGeneralUtils::MeanAndSdev(scals, mean, sdev);
		CGeneralUtils::MinMax(scals, minScal, maxScal);
		CGeneralUtils::Median(scals, 0.05, frac05);
		CGeneralUtils::Median(scals, 0.95, frac95);
		CGeneralUtils::Median(scals, 0.50, median);
		double RMS = CGeneralUtils::RMS(scals);
		
		ost << "Point scalars:" << sep;
		ost << " Min: " << spa << minScal << sep;
		ost << " Max: " << spa << maxScal << sep;
		ost << " Average: " << spa << mean << sep;
		ost << " Sdev: " << spa << sdev << sep;
		ost << " Median: " << spa << median << sep;
		ost << " RMS: " << spa << RMS << sep;
		ost << " 5% fractile: " << spa << frac05 << sep;
		ost << " 95% fractile: " << spa << frac95 << sep;
		ost << sep;
	}
	// ost << "Point data number of arrays: " << pd->GetPointData()->GetNumberOfArrays() << sep;

	scalars = vtkDoubleArray::SafeDownCast(pd->GetCellData()->GetScalars());
	if (!scalars)
	{
		if (pd->GetCellData()->GetScalars() && pd->GetCellData()->GetScalars()->GetNumberOfComponents() == 3)
		{
			ost << "Cell scalars: 3 component scalars - probably RGB" << sep;
		}
		else if (pd->GetCellData()->GetScalars())
		{
			ost << "Cell scalars: " << pd->GetCellData()->GetScalars()->GetNumberOfComponents() << " component scalars" << sep;
		}
		else
			ost << "No cell scalars" << sep;
	}
	else
	{
		std::vector<double> scals(pd->GetNumberOfPoints());

		for (int i = 0; i < pd->GetNumberOfPoints(); i++)
		{
			scals[i] = scalars->GetValue(i);
		}

		double frac05 = 0;
		double frac95 = 0;
		double median = 0;
		double mean = 0;
		double sdev = 0;
		double minScal = 0;
		double maxScal = 0;
		CGeneralUtils::MeanAndSdev(scals, mean, sdev);
		CGeneralUtils::MinMax(scals, minScal, maxScal);
		CGeneralUtils::Median(scals, 0.05, frac05);
		CGeneralUtils::Median(scals, 0.95, frac95);
		CGeneralUtils::Median(scals, 0.50, median);
		double RMS = CGeneralUtils::RMS(scals);

		ost << "Cell scalars:" << sep;
		ost << " Min: " << spa << minScal << sep;
		ost << " Max: " << spa << maxScal << sep;
		ost << " Average: " << spa << mean << sep;
		ost << " Sdev: " << spa << sdev << sep;
		ost << " Median: " << spa << median << sep;
		ost << " RMS: " << spa << RMS << sep;
		ost << " 5% fractile: " << spa << frac05 << sep;
		ost << " 95% fractile: " << spa << frac95 << sep;
		ost << sep;
	}
	//ost << "Cell data number of arrays: " << pd->GetCellData()->GetNumberOfArrays() << sep;

	vtkDataArray* norms = pd->GetPointData()->GetNormals();
	if (norms)
	{
		ost << "Contains normals" << sep;
		ost << sep;
	}
	else
	{
		ost << "No normals" << sep;
	}
	
	vtkFloatArray *TCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetTCoords());
	if (TCoords)
	{
		ost << "Contains texture coordinates" << sep;
		ost << sep;
	}
	else
	{
		ost << "No texture coordinates" << sep;
	}

	double bounds[6];
	pd->GetBounds(bounds);

	double volume = ((bounds[1]-bounds[0]) * (bounds[3]-bounds[2]) * (bounds[5]-bounds[4]));

	ost << "Bounds:" << sep
		<< " (xmin,xmax, ymin,ymax, zmin,zmax) = "
		<< "(" << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << ", " 
		<< bounds[3] << ", " << bounds[4] << ", " << bounds[5] << ")" << sep;
	ost << "Bounding box volume: " << spa << volume << sep;
	ost << sep;

	typedef std::map<int, int> CellContainer;
	CellContainer cellMap;
	for (int i = 0; i < pd->GetNumberOfCells(); i++)
	{
		cellMap[pd->GetCellType(i)]++;
	}

	ost << "Contains cell types:" << sep;
	CellContainer::const_iterator it = cellMap.begin();
	while (it != cellMap.end())
	{
		ost << "\t" << vtkCellTypes::GetClassNameFromTypeId(it->first) << " occurs "
			<< it->second << " times." << sep;
		++it;
	}

	// Now check for point data
	vtkPointData* pdata = pd->GetPointData();
	if (pdata)
	{
		ost << "Contains point data with " << pdata->GetNumberOfArrays() << " arrays."
			<< sep;
		for (int i = 0; i < pdata->GetNumberOfArrays(); i++)
		{
			ost << "\tArray " << i << " is named "
				<< (pdata->GetArrayName(i) ? pdata->GetArrayName(i) : "NULL") << sep;
		}
	}
	// Now check for cell data
	vtkCellData* cd = pd->GetCellData();
	if (cd)
	{
		ost << "Contains cell data with " << cd->GetNumberOfArrays() << " arrays."
			<< sep;
		for (int i = 0; i < cd->GetNumberOfArrays(); i++)
		{
			ost << "\tArray " << i << " is named "
				<< (cd->GetArrayName(i) ? cd->GetArrayName(i) : "NULL") << sep;
		}
	}
	// Now check for field data
	if (pd->GetFieldData())
	{
		ost << "Contains field data with " << pd->GetFieldData()->GetNumberOfArrays()
			<< " arrays." << sep;
		for (int i = 0; i < pd->GetFieldData()->GetNumberOfArrays(); i++)
		{
			if (pd->GetFieldData()->GetArray(i))
				ost << "\tArray " << i << " is named "
					<< pd->GetFieldData()->GetArray(i)->GetName() << sep;
			else
				ost << "\tArray " << i << " is NULL" << sep;
		}
	}
	ost << sep;

	return ost.str();
}

bool vtkExtMisc::CheckSurfaceConsistency( vtkPolyData *pd, std::string& msg )
{
	if (pd->GetNumberOfPoints() < 3)
	{
		msg = "Iso-surface contains less than 3 points!";
		return false;
	}

	vtkFeatureEdges *features = vtkFeatureEdges::New();
	features->SetInputData(pd);
	features->SetNonManifoldEdges(1);
	features->SetBoundaryEdges(0);
	features->SetFeatureEdges(0);
	features->SetManifoldEdges(0);
	features->Update();

	int nfeat = features->GetOutput()->GetNumberOfLines();
	features->Delete();

	if (nfeat > 0)
	{
		std::ostringstream ost;
		ost << "Surface contains " << nfeat << " non-manifold edges!";
		msg = ost.str();
		return false;
	}
	return true;
}

vtkPolyData * vtkExtMisc::MultiReadSurface( const std::string &fname )
{
	vtkPolyData *pdout = vtkPolyData::New();

	std::string ext = CGeneralUtils::GetExtensionFromFilename(fname);

	if (ext == "vtk")
	{
		vtkPolyDataReader *reader = vtkExtMisc::SafeReadPolyData(fname);
		if (!reader)
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "vtp")
	{
		vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (!reader)
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "stl")
	{
		vtkSTLReader *reader = vtkExtMisc::SafeReadPolyDataSTL(fname);
		if (!reader)
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "obj")
	{
		vtkOBJReader *reader = vtkExtMisc::SafeReadPolyDataOBJ(fname);
		if (!reader)
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	//else if (ext == "aranz")
	//{
	//	vtkARANZReader *reader = vtkARANZReader::New();
	//	reader->SetFileName(fname.c_str());
	//	reader->Update();
	//	if (!reader || reader->GetOutput()->GetNumberOfPoints() < 1)
	//	{
	//		pdout->Delete();
	//		return NULL;
	//	}
	//	pdout->DeepCopy(reader->GetOutput());
	//	reader->Delete();
	//	return pdout;
	//}
	else if (ext == "ply")
	{
		vtkPLYReader *reader = vtkPLYReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (!reader || reader->GetOutput()->GetNumberOfPoints() < 1)
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	//else if (ext == "3dmtxt")
	//{
	//	vtk3DMDTxtReader *reader = vtk3DMDTxtReader::New();
	//	reader->SetFileName(fname.c_str());
	//	reader->Update();
	//	if (!reader || reader->GetOutput()->GetNumberOfPoints() < 1)
	//	{
	//		pdout->Delete();
	//		return NULL;
	//	}
	//	pdout->DeepCopy(reader->GetOutput());
	//	reader->Delete();
	//	return pdout;
	//}
	else if (ext == "dat")
	{
		vtkPolyDataTextReader *reader = vtkPolyDataTextReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (reader->ReaderError())
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "csv")
	{
		vtkPolyDataTextReader *reader = vtkPolyDataTextReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (reader->ReaderError())
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "pts")
	{
		vtkPolyDataTextReader *reader = vtkPolyDataTextReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (reader->ReaderError())
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "txt")
	{
		vtkPolyDataTextReader *reader = vtkPolyDataTextReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (reader->ReaderError())
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "bnd" || ext == "pse") // BU-3DFE landmarks
	{
		vtkPolyDataTextReader *reader = vtkPolyDataTextReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (reader->ReaderError())
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "anno") // Jens Fagertun simple landmarks
	{
		vtkPolyDataTextReader *reader = vtkPolyDataTextReader::New();
		reader->SetFileName(fname.c_str());
		reader->Update();
		if (reader->ReaderError())
		{
			pdout->Delete();
			return NULL;
		}
		pdout->DeepCopy(reader->GetOutput());
		reader->Delete();
		return pdout;
	}
	else if (ext == "pp") // MeshLab XML Landmark file
	{
		CLandmarkCollection LMCollection;
		if (!LMCollection.ReadFromXLMFile(fname))
		{
			pdout->Delete();
			return NULL;
		}
		vtkPolyData *pd = LMCollection.GetAsVTKPolydata();
		pdout->DeepCopy(pd);
		pd->Delete();
		return pdout;
	}
	std::cout << "Unknown extension " << ext << std::endl;

	// No good reader found
	pdout->Delete();
	return NULL;
}


bool vtkExtMisc::MultiWriteSurface(const std::string& fname, vtkPolyData* pd, bool WriteNormals, bool WriteScalars)
{
	std::string ext = CGeneralUtils::GetExtensionFromFilename(fname);

	bool IsScalars = false;
	bool IsRBGPerVertex = false;

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetScalars());
	if (!scalars)
	{
		IsScalars = true;
		if (pd->GetPointData()->GetScalars() && pd->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
		{
			IsRBGPerVertex = true;
		}
	}
	
	if (ext == "vtk")
	{
		return vtkExtMisc::WritePDVTK(pd, fname);
	}
	else if (ext == "vtp")
	{
		vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
		writer->SetInputData(pd);
		writer->SetFileName(fname.c_str());
		bool result = (writer->Write() == 1);
		writer->Delete();
		return result;
	}
	else if (ext == "stl")
	{
		return vtkExtMisc::WritePDSTL(pd, fname);
	}
	//else if (ext == "aranz")
	//{
	//	vtkARANZWriter *writer = vtkARANZWriter::New();
	//	writer->SetInputData(pd);
	//	writer->SetFileName(fname.c_str());
	//	bool result = (writer->Write() == 1);
	//	writer->Delete();
	//	return result;
	//}
	else if (ext == "ply")
	{
		vtkPLYWriter *writer = vtkPLYWriter::New();
		writer->SetInputData(pd);
		if (IsRBGPerVertex)
		{
			writer->SetArrayName("Colors");
			writer->SetColorModeToDefault();
		}
		writer->SetFileName(fname.c_str());
		bool result = (writer->Write() == 1);
		writer->Delete();
		return result;
	}
	else if (ext == "obj")
	{
		vtkOBJWriter *writer = vtkOBJWriter::New();
		writer->SetInputData(pd);
		writer->SetFileName(fname.c_str());
		writer->Update();
		return true;
	}
	else if (ext == "txt")
	{
		vtkRawPolyWriter *writer = vtkRawPolyWriter::New();
		writer->SetInputData(pd);
		writer->SetFileName(fname.c_str());
		writer->SetWriteNormals(WriteNormals);
		writer->SetWriteScalars(WriteScalars);
		bool result = (writer->Write() == 1);
		writer->Delete();
		return result;
	}
	std::cout << "Unknown extension " << ext << std::endl;
	return false;
}
