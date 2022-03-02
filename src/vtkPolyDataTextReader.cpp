#include "vtkPolyDataTextReader.h"

#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkMath.h"

#include "GeneralUtils.h"
#include <fstream>
#include <sstream>

vtkStandardNewMacro(vtkPolyDataTextReader);

// Description:
// Instantiate object with NULL filename.
vtkPolyDataTextReader::vtkPolyDataTextReader()
{
	this->FileName = NULL;
	
	this->SetNumberOfInputPorts(0);
	ErrorReading = false;
}

vtkPolyDataTextReader::~vtkPolyDataTextReader()
{
	if (this->FileName)
	{
		delete [] this->FileName;
		this->FileName = NULL;
	}
}

int vtkPolyDataTextReader::RequestData(
								vtkInformation *vtkNotUsed(request),
								vtkInformationVector **vtkNotUsed(inputVector),
								vtkInformationVector *outputVector)
{
	ErrorReading = false;
//	int i;

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	
	// get the ouptut
	vtkPolyData *output = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));
	
	if (!this->FileName) 
	{
		vtkErrorMacro(<< "A FileName must be specified.");
		ErrorReading = true;
		return 0;
	}

	bool bndMode = false; // Are we in the mode of "BU-3D" data with a vertex id first
	bool psemode = false;
	if (CGeneralUtils::GetExtensionFromFilename(std::string(this->FileName)) == "bnd")
	{
		bndMode = true;
	}
	else if (CGeneralUtils::GetExtensionFromFilename(std::string(this->FileName)) == "pse")  // Last entry is a normal
	{
		psemode = true;
	}

	// Initialize
	std::ifstream istr(this->FileName);
	if (!istr)
	{
		vtkErrorMacro(<< "File " << this->FileName << " not found");
		ErrorReading = true;
		return 0;
	}

	// Find out if there are delimeters and how many number per line
	// 3 = points only
	// 4 = points + scalars
	// 6 = points + normals
	// 7 = points + normals + scalars
	// 9 = points + normals + RGB (you can use 0 0 0 as dummy normals)
	std::string str;
	std::getline(istr, str);
	std::string delim = " ";
	bool SpaceDelim = true;

	size_t pos = str.find(",");
	if (pos != -1)
	{
		delim = ",";
		SpaceDelim = false;
	}
	else
	{
		size_t pos = str.find(";");
		if (pos != -1)
		{
			delim = ";";
			SpaceDelim = false;
		}
		else
		{
			size_t pos = str.find("\t");
			if (pos != -1)
			{
				delim = "\t";
				SpaceDelim = false;
			}
		}
	}

	std::vector<double> nums;
	int nsp = 0;
	if (SpaceDelim)
		nsp = CGeneralUtils::GetNumbersFromString(str, nums);
	else
		nsp = CGeneralUtils::SplitStringIntoNumbers(str, delim, nums);

	if (nsp < 2 || nsp > 8)
	{
		ErrorReading = true;
		return 0;
	}

	// X,Y,Z
	if (nsp == 2)
	{
		vtkPoints *pts = vtkPoints::New();
		vtkCellArray *verts = vtkCellArray::New();

		bool stop = false;
		do
		{
			// We already read the first line
			std::vector<double> st;
			int nst = 0;
			if (SpaceDelim)
				nst = CGeneralUtils::GetNumbersFromString(str, st);
			else
				nst = CGeneralUtils::SplitStringIntoNumbers(str, delim, st);

			if (nst != nsp)
			{
				stop = true;
			}
			else
			{
				double p[3];
				p[0] = st[0];
				p[1] = st[1];
				p[2] = st[2];

				vtkIdType id = pts->InsertNextPoint(p);
				verts->InsertNextCell(1);
				verts->InsertCellPoint(id);

				// get next line
				std::getline(istr, str);

				if (istr.fail())
				{
					stop = true;
				}

			}
		}
		while (!stop);

		output->SetPoints(pts);
		output->SetVerts(verts);
		pts->Delete();
		verts->Delete();
		if (output->GetNumberOfPoints() < 1)
		{
			output->Delete();
			ErrorReading = true;
			return 0;
		}
	}

	// X,Y,Z, scalars
	// or if bndMode = true (vertexID, X, Y, Z)
	if (nsp == 3)
	{
		vtkPoints *pts = vtkPoints::New();
		vtkCellArray *verts = vtkCellArray::New();
		vtkDoubleArray *scalars = vtkDoubleArray::New();
		scalars->SetNumberOfComponents(1);

		bool stop = false;
		do
		{
			if (psemode)
			{
				size_t res = str.find("n", 0);
				if (res >= 0 )
					stop = true;
			}

			if (!stop)
			{
				// We already read the first line
				std::vector<double> st;
				int nst = 0;
				if (SpaceDelim)
					nst = CGeneralUtils::GetNumbersFromString(str, st);
				else
					nst = CGeneralUtils::SplitStringIntoNumbers(str, delim, st);

				if (nst != nsp)
				{
					stop = true;
				}
				else
				{
					double p[3];
					if (bndMode || psemode)
					{
						p[0] = st[1];
						p[1] = st[2];
						p[2] = st[3];
					}
					else
					{
						p[0] = st[0];
						p[1] = st[1];
						p[2] = st[2];
					}

					vtkIdType id = pts->InsertNextPoint(p);
					verts->InsertNextCell(1);
					verts->InsertCellPoint(id);

					if (!bndMode && !psemode)
					{
						double scal = st[3];
						scalars->InsertNextValue(scal);
					}

					// get next line
					std::getline(istr, str);

					if (istr.fail())
					{
						stop = true;
					}
				}
			}
		}
		while (!stop);

		output->SetPoints(pts);
		output->SetVerts(verts);
		if (!bndMode && !psemode)
		{
			output->GetPointData()->SetScalars(scalars);
		}

		pts->Delete();
		verts->Delete();
		scalars->Delete();
		if (output->GetNumberOfPoints() < 1)
		{
			output->Delete();
			ErrorReading = true;
			return 0;
		}
	}

	// X,Y,Z, Xn, Yn, Zn  (normal)
	if (nsp == 5)
	{
		vtkPoints *pts = vtkPoints::New();
		vtkCellArray *verts = vtkCellArray::New();
		vtkDoubleArray *Normals = vtkDoubleArray::New();
		Normals->SetNumberOfComponents(3);
		Normals->SetName("Normals");

		bool stop = false;
		do
		{
			// We already read the first line
			std::vector<double> st;
			int nst = 0;
			if (SpaceDelim)
				nst = CGeneralUtils::GetNumbersFromString(str, st);
			else
				nst = CGeneralUtils::SplitStringIntoNumbers(str, delim, st);

			if (nst != nsp)
			{
				stop = true;
			}
			else
			{
				double p[3];
				p[0] = st[0];
				p[1] = st[1];
				p[2] = st[2];

				vtkIdType id = pts->InsertNextPoint(p);
				verts->InsertNextCell(1);
				verts->InsertCellPoint(id);

				double n[3];
				n[0] = st[3];
				n[1] = st[4];
				n[2] = st[5];

				vtkMath::Normalize(n);
				Normals->InsertNextTuple(n);

				// get next line
				std::getline(istr, str);

				if (istr.fail())
				{
					stop = true;
				}

			}
		}
		while (!stop);

		output->SetPoints(pts);
		output->SetVerts(verts);
		output->GetPointData()->SetNormals(Normals);

		Normals->Delete();
		pts->Delete();
		verts->Delete();
		if (output->GetNumberOfPoints() < 1)
		{
			output->Delete();
			ErrorReading = true;
			return 0;
		}
	}

	// X,Y,Z, Xn, Yn, Zn  (normal), scalar
	if (nsp == 6)
	{
		vtkPoints *pts = vtkPoints::New();
		vtkCellArray *verts = vtkCellArray::New();
		vtkDoubleArray *Normals = vtkDoubleArray::New();
		Normals->SetNumberOfComponents(3);
		Normals->SetName("Normals");
		vtkDoubleArray *scalars = vtkDoubleArray::New();
		scalars->SetNumberOfComponents(1);

		bool stop = false;
		do
		{
			// We already read the first line
			std::vector<double> st;
			int nst = 0;
			if (SpaceDelim)
				nst = CGeneralUtils::GetNumbersFromString(str, st);
			else
				nst = CGeneralUtils::SplitStringIntoNumbers(str, delim, st);

			if (nst != nsp)
			{
				stop = true;
			}
			else
			{
				double p[3];
				p[0] = st[0];
				p[1] = st[1];
				p[2] = st[2];

				vtkIdType id = pts->InsertNextPoint(p);
				verts->InsertNextCell(1);
				verts->InsertCellPoint(id);

				double n[3];
				n[0] = st[3];
				n[1] = st[4];
				n[2] = st[5];

				vtkMath::Normalize(n);
				Normals->InsertNextTuple(n);

				double scal = st[6];
				scalars->InsertNextValue(scal);

				// get next line
				std::getline(istr, str);

				if (istr.fail())
				{
					stop = true;
				}

			}
		}
		while (!stop);

		output->SetPoints(pts);
		output->SetVerts(verts);
		output->GetPointData()->SetNormals(Normals);
		output->GetPointData()->SetScalars(scalars);

		scalars->Delete();
		Normals->Delete();
		pts->Delete();
		verts->Delete();
		if (output->GetNumberOfPoints() < 1)
		{
			output->Delete();
			ErrorReading = true;
			return 0;
		}
	}

	// X,Y,Z, Xn, Yn, Zn  (normal), R, G, B
	if (nsp == 8)
	{
		vtkPoints *pts = vtkPoints::New();
		vtkCellArray *verts = vtkCellArray::New();
		vtkDoubleArray *Normals = vtkDoubleArray::New();
		Normals->SetNumberOfComponents(3);
		Normals->SetName("Normals");

		vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("Colors");

		//vtkDoubleArray *colors = vtkDoubleArray::New();
		//colors->SetNumberOfComponents(3);
		//colors->SetName("Colors");

		bool reallyNormals = true;
		bool stop = false;
		do
		{
			// We already read the first line
			std::vector<double> st;
			int nst = 0;
			if (SpaceDelim)
				nst = CGeneralUtils::GetNumbersFromString(str, st);
			else
				nst = CGeneralUtils::SplitStringIntoNumbers(str, delim, st);


			if (nst != nsp)
			{
				stop = true;
			}
			else
			{
				double p[3];
				p[0] = st[0];
				p[1] = st[1];
				p[2] = st[2];

				vtkIdType id = pts->InsertNextPoint(p);
				verts->InsertNextCell(1);
				verts->InsertCellPoint(id);

				double n[3];
				n[0] = st[3];
				n[1] = st[4];
				n[2] = st[5];

				// Check for dummy normals
				if (reallyNormals && n[0] == 0 && n[1] == 0 && n[2] == 0)
				{
					reallyNormals = false;
				}
				if (reallyNormals)
				{
					vtkMath::Normalize(n);
					Normals->InsertNextTuple(n);
				}

				unsigned char RGB[3];

//				double RGB[3];
				RGB[0] = st[6] * 255;
				RGB[1] = st[7] * 255;
				RGB[2] = st[8] * 255;
				colors->InsertNextTypedTuple(RGB);

				// get next line
				std::getline(istr, str);

				if (istr.fail())
				{
					stop = true;
				}

			}
		}
		while (!stop);

		output->SetPoints(pts);
		output->SetVerts(verts);

		if (reallyNormals)
		{
			output->GetPointData()->SetNormals(Normals);
		}

		output->GetPointData()->SetScalars(colors);
		
		colors->Delete();
		Normals->Delete();
		pts->Delete();
		verts->Delete();
		if (output->GetNumberOfPoints() < 1)
		{
			output->Delete();
			ErrorReading = true;
			return 0;
		}
	}

	return 1;
}

void vtkPolyDataTextReader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
	
	os << indent << "File Name: " 
		<< (this->FileName ? this->FileName : "(none)") << "\n";
	
}

bool vtkPolyDataTextReader::ReaderError() const
{
	return ErrorReading;
}
