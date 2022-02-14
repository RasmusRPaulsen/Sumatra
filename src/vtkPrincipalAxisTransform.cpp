/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPrincipalAxisTransform.cpp,v $
  Language:  C++
  Date:      $Date: 2003/12/29 16:13:21 $
  Version:   $Revision: 1.8 $
  
  Made by Rasmus Paulsen
  email:  rrp@imm.dtu.dk
  web:    www.imm.dtu.dk/~rrp/VTK

  This class is not mature enough to enter the official VTK release.
=========================================================================*/
#include "vtkPrincipalAxisTransform.h"
#include "vtkPoints.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkPrincipalAxisTransform);

//----------------------------------------------------------------------------

vtkPrincipalAxisTransform::vtkPrincipalAxisTransform()
: vtkTransform()
{
	this->Source = NULL;
	this->DoRotate = 1;
	this->DoTranslate= 1;
	this->DoScale = 1;
}

//----------------------------------------------------------------------------

vtkPrincipalAxisTransform::~vtkPrincipalAxisTransform()
{
	ReleaseSource();
}

//----------------------------------------------------------------------------

void vtkPrincipalAxisTransform::SetSource(vtkPoints *source)
{
	if (this->Source == source)
	{
		return;
	}
	
	if (this->Source)
	{
		this->ReleaseSource();
	}
	
	if (source)
	{
		source->Register(this);
	}
	
	this->Source = source;
	this->Modified();
}

//----------------------------------------------------------------------------

void vtkPrincipalAxisTransform::ReleaseSource(void) {
	if (this->Source) 
	{
		this->Source->UnRegister(this);
		this->Source = NULL;
	}
}

//------------------------------------------------------------------------

vtkMTimeType vtkPrincipalAxisTransform::GetMTime()
{
	unsigned long result = this->vtkLinearTransform::GetMTime();
	unsigned long mtime;
	
	if (this->Source)
	{
		mtime = this->Source->GetMTime(); 
		if (mtime > result)
		{
			result = mtime;
		}
	}
	
	return result;
}


//----------------------------------------------------------------------------

//void vtkPrincipalAxisTransform::Inverse()
//{
//}


//----------------------------------------------------------------------------

vtkAbstractTransform *vtkPrincipalAxisTransform::MakeTransform()
{
	return vtkPrincipalAxisTransform::New(); 
}

//----------------------------------------------------------------------------

void vtkPrincipalAxisTransform::InternalDeepCopy(vtkAbstractTransform *transform)
{
	vtkPrincipalAxisTransform *t = (vtkPrincipalAxisTransform *)transform;
	
	this->SetSource(t->GetSource());
	this->Modified();
}

//----------------------------------------------------------------------------

void vtkPrincipalAxisTransform::InternalUpdate()
{
	int i,j;

	vtkIdType pointId;

	// Check source, target
	
	if (this->Source == NULL || !this->Source->GetNumberOfPoints())
	{
		vtkErrorMacro(<<"Can't execute with NULL or empty input");
		return;
	}
	
	//
	// Compute mean
	//
	vtkIdType numPts = Source->GetNumberOfPoints();

	double mean[3] = {0.0, 0.0, 0.0};
	for (pointId = 0; pointId < numPts; pointId++)
	{
		double *x = Source->GetPoint(pointId);
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
		double *x = Source->GetPoint(pointId);

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
	
	Identity();
	if (DoRotate)
	{
		if (DoScale)
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					this->Matrix->SetElement(i, j, v[i][j] * sqrt(evals[j]));
				}
			}
		}
		else
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					this->Matrix->SetElement(i, j, v[i][j]);
				}
			}
		}
	}
	else if (DoScale)
	{
		for (i = 0; i < 3; i++)
		{
			this->Matrix->SetElement(i, i, sqrt(evals[i]));
		}
	}

	if (DoTranslate)
	{
	//	Translate(mean[0], mean[1], mean[2]);

		// Hard code translation
		this->Matrix->SetElement(0, 3, mean[0]);
		this->Matrix->SetElement(1, 3, mean[1]);
		this->Matrix->SetElement(2, 3, mean[2]);
	}

	// Check
	double (*matrix)[4] = this->Matrix->Element;
	double U[3][3], VT[3][3];
	
	for (i = 0; i < 3; i++) 
	{
		U[0][i] = matrix[0][i];
		U[1][i] = matrix[1][i];
		U[2][i] = matrix[2][i];
	}
	
	double scale[3];
	vtkMath::SingularValueDecomposition3x3(U, U, scale, VT);

//	cout << "check scale " << scale[0] << " " << scale[1] << " " << scale[2] << endl;

	// This is dirty, but I do not know how else to avoid the scale being negative!!!
	if (scale[0] < 0 || scale[1] < 0 || scale[1] < 0)
	{
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				this->Matrix->SetElement(i, j, -(this->Matrix->GetElement(i,j)));
			}
		}
	}
}

//----------------------------------------------------------------------------

void vtkPrincipalAxisTransform::PrintSelf(ostream& os, vtkIndent indent)
{
	vtkLinearTransform::PrintSelf(os,indent);
	
	if ( this->Source ) 
	{
		os << indent << "Source: " << this->Source << "\n";
	}
	else 
	{
		os << indent << "Source: (none)\n";
	}
}
