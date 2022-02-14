#include "vtkRawPolyWriter.h"

#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkByteSwap.h"
#include "vtkCellArray.h"
#include "vtkErrorCode.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"

#if !defined(_WIN32) || defined(__CYGWIN__)
# include <unistd.h> /* unlink */
#else
# include <io.h> /* unlink */
#endif

vtkStandardNewMacro(vtkRawPolyWriter);

vtkRawPolyWriter::vtkRawPolyWriter()
{
	this->FileType = VTK_ASCII;
	this->WriteScalars = 0;
	this->WriteNormals = 0;
}

void vtkRawPolyWriter::WriteData()
{
	vtkPoints *pts;
	vtkCellArray *polys;
	vtkPolyData *input = this->GetInput();
	
	polys = input->GetPolys();
	pts = input->GetPoints();
	if (pts == NULL || polys == NULL )
	{
		vtkErrorMacro(<<"No data to write!");
		return;
	}
	
	if ( this->FileName == NULL)
	{
		vtkErrorMacro(<< "Please specify FileName to write");
		this->SetErrorCode(vtkErrorCode::NoFileNameError);
		return;
	}
	
	{
		this->WriteAsciiRAW(input);
		if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
		{
			vtkErrorMacro("Ran out of disk space; deleting file: "
				<< this->FileName);
			_unlink(this->FileName);
		}
	}
}

static char header[]="Visualization Toolkit generated RAW File                                        ";

void vtkRawPolyWriter::WriteAsciiRAW(vtkPolyData *pd)
{
	vtkCellArray *polys = pd->GetPolys();
	vtkPoints *pts = pd->GetPoints();

	FILE *fp;
	
	if ((fp = fopen(this->FileName, "w")) == NULL)
	{
		vtkErrorMacro(<< "Couldn't open file: " << this->FileName);
		this->SetErrorCode(vtkErrorCode::CannotOpenFileError);
		return;
	}
	
	const int NPoints = pts->GetNumberOfPoints(); 

	vtkDoubleArray *scalars = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetScalars());
	bool IsScalars = (scalars != NULL);

	vtkDataArray *normals = pd->GetPointData()->GetNormals(); 
	bool isNormals = (normals != NULL);

	for (int i = 0; i < NPoints; i++)
	{
		double v[3];
		pts->GetPoint(i, v);

		if (fprintf (fp, "%.8g %.8g %.8g", v[0], v[1], v[2]) < 0)
		{
			fclose(fp);
			this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
			return;
		}
		if (isNormals && WriteNormals)
		{
			double n[3];
			normals->GetTuple(i, n);
			if (fprintf (fp, " %.8g %.8g %.8g", n[0], n[1], n[2]) < 0)
			{
				fclose(fp);
				this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
				return;
			}
		}
		if (IsScalars && WriteScalars)
		{
			double s = scalars->GetValue(i);
			if (fprintf (fp, " %.8g", s) < 0)
			{
				fclose(fp);
				this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
				return;
			}
		}
		if (fprintf (fp, "\n") < 0)
		{
			fclose(fp);
			this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
			return;
		}
	}


	//if (IsScalars)
	//{
	//	for (int i = 0; i < NPoints; i++)
	//	{
	//		double s = scalars->GetValue(i);
	//		double v[3];
	//		pts->GetPoint(i, v);

	//		if (fprintf (fp, "%.8g %.8g %.8g %.8g\n", v[0], v[1], v[2],s) < 0)
	//		{
	//			fclose(fp);
	//			this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
	//			return;
	//		}
	//	}
	//}
	//else
	//{
	//	for (int i = 0; i < NPoints; i++)
	//	{
	//		double v[3];
	//		pts->GetPoint(i, v);

	//		if (fprintf (fp, "%.8g %.8g %.8g\n", v[0], v[1], v[2]) < 0)
	//		{
	//			fclose(fp);
	//			this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
	//			return;
	//		}
	//	}
	//}

	fclose (fp);
}

//----------------------------------------------------------------------------
void vtkRawPolyWriter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
