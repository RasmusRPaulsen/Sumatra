#include "SurfaceProperties.h"

#include <sys/stat.h> // for file exists

#include <vtkJPEGReader.h>
#include <vtkTIFFReader.h>
#include <vtkTexture.h>
#include <vtkBMPReader.h>
#include <vtkPNGReader.h>
#include <vtkMatrix4x4.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLookupTable.h>
//#include <vtkARANZWriter.h>
//#include <vtkARANZReader.h>
#include <vtkRawPolyWriter.h>
#include <vtkPLYReader.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
//#include "vtk3DMDTxtReader.h"

#include "vtkExtMisc.h"
#include "GeneralUtils.h"
#include "vtkVRMLImporter.h"

 CSurfaceProperties::CSurfaceProperties()
{
}

CSurfaceProperties::CSurfaceProperties(vtkLookupTable *lut, CSumatraSettings* settings)
{
	m_shortname = "";
	m_fullname = "";
	m_UndoAvailable = false;
	m_UndoPolyData = NULL;
	m_UndoMatrix = NULL;
	m_polyData = NULL;
	m_Texture = NULL;
	m_mapper = NULL;
	m_actor = NULL;
	m_lookup = lut;
	mSettings = settings;
}


CSurfaceProperties::~CSurfaceProperties()
{
	if (m_polyData)
		m_polyData->Delete();
	if (m_mapper)
		m_mapper->Delete();
	if (m_Texture)
		m_Texture->Delete();
	if (m_actor)
		m_actor->Delete();
	if (m_UndoMatrix)
		m_UndoMatrix->Delete();
	if (m_UndoPolyData)
		m_UndoPolyData->Delete();
}

void CSurfaceProperties::InitialiseSurface(vtkPolyData *pd)
{
/*
Ear:  255 210 150
Probe: 200 200 200 
Rays: 0 70 255
Points 0 0 100

*/
	m_polyData = vtkPolyData::New();
	m_polyData->DeepCopy(pd);

	m_mapper = vtkPolyDataMapper::New();
	m_mapper->SetInputData(m_polyData);

	m_actor = vtkActor::New();
	m_actor->SetMapper(m_mapper);

	double color[3];
	// m_colorManager->GetNextColor(color);
	mSettings->GetNextColor(color);

	vtkProperty *prop = m_actor->GetProperty();
//	prop->SetColor(255.0/255.0, 210.0/255.0, 150.0/255.0);
	prop->SetColor(color);
	prop->SetSpecularColor(1, 1, 1);
	prop->SetSpecular(0.3);
	prop->SetSpecularPower(20);
	prop->SetAmbient(0.2);
	prop->SetDiffuse(0.8);
	prop->SetOpacity(1.0);
	prop->SetPointSize(mSettings->mPointSize);
	prop->SetLineWidth(mSettings->mLineWidth);
	prop->SetEdgeColor(0, 0, 0);
	// prop->SetVertexVisibility(true);

	UpdateScalarProperties();

	// Check for RGB scalar values
	if (m_polyData->GetPointData()->GetScalars() && m_polyData->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
	{
		prop->SetColor(1,1,1);
		prop->SetSpecularColor(1, 1, 1);
		prop->SetSpecular(0);
		prop->SetSpecularPower(0);
		prop->SetAmbient(1.0);
		prop->SetDiffuse(0);
	}
}

void CSurfaceProperties::UpdateScalarProperties()
{
	// no scalars -> do not update lookup table
	if (m_lookup != NULL && m_polyData->GetPointData()->GetScalars() != NULL)
	{
		double range[2];
		m_polyData->GetPointData()->GetScalars()->GetRange(range);

		if (range[0] != range[1])
		{
			m_mapper->SetScalarRange(range);
			m_mapper->SetLookupTable(m_lookup);
			m_mapper->SetScalarVisibility(true);
		}
		else
		{
			m_mapper->SetScalarVisibility(false);
		}
	}
	else
	{
		m_mapper->SetScalarVisibility(false);
	}
}

void CSurfaceProperties::StoreState()
{
	if (m_UndoPolyData == NULL)
		m_UndoPolyData = vtkPolyData::New();
	if (m_UndoMatrix == NULL)
		m_UndoMatrix = vtkMatrix4x4::New();

	m_UndoPolyData->DeepCopy(m_polyData);
	m_UndoMatrix->DeepCopy(m_actor->GetMatrix());
	m_UndoName1 = m_shortname;
	m_UndoName2 = m_fullname;
	m_UndoAvailable = true;
}

void CSurfaceProperties::RestoreState()
{
	if (m_UndoPolyData == NULL)
		return;

	m_polyData->DeepCopy(m_UndoPolyData);
	m_actor->GetMatrix()->DeepCopy(m_UndoMatrix);
	m_shortname = m_UndoName1;
	m_fullname  = m_UndoName2;
	m_UndoAvailable = false;
}

bool CSurfaceProperties::UndoAvailable() const
{
	return m_UndoAvailable;
}

bool CSurfaceProperties::HaveScalars() const
{
	if (m_polyData->GetPointData()->GetScalars() == NULL)
		return false;
	return true;
}

bool CSurfaceProperties::GetScalarRange(double *range)
{
	// no scalars -> do not update lookup table
	if (m_polyData->GetPointData()->GetScalars() == NULL)
		return false;
		
//	m_polyData->GetPointData()->GetScalars()->GetRange(range);
	m_mapper->GetScalarRange(range);

	if (range[0] != range[1])
		return false;

	return true;
}

void CSurfaceProperties::SetScalarRange(double *range)
{
	m_mapper->SetScalarRange(range);
	m_mapper->Update();
}

bool CSurfaceProperties::SaveToFile(const std::string& fname, bool ApplyUserTransform, bool Binary, bool WriteNormals, bool WriteScalars)
{
	std::string ext = CGeneralUtils::GetExtensionFromFilename(fname);

	m_fullname = fname;
	if (ApplyUserTransform)
	{
		vtkTransform *pdTrans = vtkTransform::New();
		 pdTrans->SetMatrix(m_actor->GetMatrix());

		vtkTransformPolyDataFilter *tf= vtkTransformPolyDataFilter::New();
		 tf->SetInputData(m_polyData);
		 tf->SetTransform(pdTrans);
		 tf->Update();

		bool result = false;
		if (ext == "vtk")
		{
			 result = vtkExtMisc::WritePDVTK(tf->GetOutput(), fname, Binary);
		}
		else if (ext == "stl")
		{
			result = vtkExtMisc::WritePDSTL(tf->GetOutput(), fname, Binary);
		}
		//else if (ext == "aranz")
		//{
		//	vtkARANZWriter *writer = vtkARANZWriter::New();
		//	 writer->SetInputConnection(tf->GetOutputPort());
		//	 writer->SetFileName(fname.c_str());
		//	 result = (writer->Write() == 1);
		//	 writer->Delete();
		//}
		else if (ext == "txt")
		{
			vtkRawPolyWriter *writer = vtkRawPolyWriter::New();
			writer->SetInputConnection(tf->GetOutputPort());
			writer->SetFileName(fname.c_str());
			writer->SetWriteNormals(WriteNormals);
			writer->SetWriteScalars(WriteScalars);
			result = (writer->Write() == 1);
			writer->Delete();
		}
		else
		{
			result = vtkExtMisc::MultiWriteSurface(fname, tf->GetOutput(), WriteNormals);
		}
		tf->Delete();
		pdTrans->Delete();
		return result;
	}
	else
	{
		if (ext == "vtk")
		{
			return vtkExtMisc::WritePDVTK(m_polyData, fname, Binary);
		}
		else if (ext == "stl")
		{
			return vtkExtMisc::WritePDSTL(m_polyData, fname, Binary);
		}
		//else if (ext == "aranz")
		//{
		//	vtkARANZWriter *writer = vtkARANZWriter::New();
		//	 writer->SetInputData(m_polyData);
		//	 writer->SetFileName(fname.c_str());
		//	 bool result = (writer->Write() == 1);
		//	 writer->Delete();
		//	 return result;
		//}
		else
		{
			return vtkExtMisc::MultiWriteSurface(fname, m_polyData, WriteNormals);
		}
	}
	return false;
}


void CSurfaceProperties::SetTexture( vtkImageData *tex )
{
	m_Texture = vtkTexture::New();
	m_Texture->SetInterpolate(1);
	m_Texture->SetInputData((vtkDataObject*)tex);

	m_actor->SetTexture(m_Texture);

	// Reset color to pure white
	m_actor->GetProperty()->SetColor(1,1,1);
	m_actor->GetProperty()->SetAmbient(1.0);
}


bool CSurfaceProperties::ReadTexture(const std::string& fname)
{
	// First check for texture coordinates in the polydata
	vtkDataArray *TexCoords = m_polyData->GetPointData()->GetTCoords();
	if (!TexCoords)
		return false;

	std::string PNGName = CGeneralUtils::StripExtensionFromFilename(fname) + ".png";
	std::string TIFName = CGeneralUtils::StripExtensionFromFilename(fname) + ".tif";
	std::string BMPName = CGeneralUtils::StripExtensionFromFilename(fname) + ".bmp";
	std::string JPGName = CGeneralUtils::StripExtensionFromFilename(fname) + ".jpg";


	if (file_exists(PNGName.c_str()))
	{
		vtkPNGReader *textureImage = vtkPNGReader::New();
		textureImage->SetFileName(PNGName.c_str());
		textureImage->Update();

		m_Texture = vtkTexture::New();
		m_Texture->SetInterpolate(1);
		m_Texture->SetInputConnection(textureImage->GetOutputPort());

		textureImage->Delete();

		m_actor->SetTexture(m_Texture);
	
		// Reset color to pure white
		m_actor->GetProperty()->SetColor(1,1,1);
		m_actor->GetProperty()->SetAmbient(1.0);

		return true;
	}
	if (file_exists(TIFName.c_str()))
	{
		vtkTIFFReader *textureImage = vtkTIFFReader::New();
		textureImage->SetFileName(TIFName.c_str());
		textureImage->Update();

		m_Texture = vtkTexture::New();
		m_Texture->SetInterpolate(1);
		m_Texture->SetInputConnection(textureImage->GetOutputPort());

		textureImage->Delete();

		m_actor->SetTexture(m_Texture);

		// Reset color to pure white
		m_actor->GetProperty()->SetColor(1,1,1);
		m_actor->GetProperty()->SetAmbient(1.0);

		return true;
	}

	if (file_exists(BMPName.c_str()))
	{
		vtkBMPReader *textureImage = vtkBMPReader::New();
		textureImage->SetFileName(BMPName.c_str());
		textureImage->Update();

		m_Texture = vtkTexture::New();
		m_Texture->SetInterpolate(1);
		m_Texture->SetInputConnection(textureImage->GetOutputPort());

		textureImage->Delete();

		m_actor->SetTexture(m_Texture);

		// Reset color to pure white
		m_actor->GetProperty()->SetColor(1,1,1);
		m_actor->GetProperty()->SetAmbient(1.0);

		return true;
	}
	if (file_exists(JPGName.c_str()))
	{
		vtkJPEGReader *textureImage = vtkJPEGReader::New();
		textureImage->SetFileName(JPGName.c_str());
		textureImage->Update();

		m_Texture = vtkTexture::New();
		m_Texture->SetInterpolate(1);
		//m_Texture->SetQualityTo32Bit();
		//m_Texture->SetRestrictPowerOf2ImageSmaller(0);
		m_Texture->SetInputConnection(textureImage->GetOutputPort());

		textureImage->Delete();

		m_actor->SetTexture(m_Texture);

		// Reset color to pure white
		m_actor->GetProperty()->SetColor(1,1,1);
		m_actor->GetProperty()->SetAmbient(1.0);

		return true;
	}

	return true;
}

bool CSurfaceProperties::ReadFromFile(const std::string& fname)
{
	m_fullname = fname;
	m_shortname = CGeneralUtils::StripPathAndExtensionFromFilename(fname);

	std::string ext = CGeneralUtils::GetExtensionFromFilename(fname);

	if (ext == "vtk")
	{
		vtkPolyDataReader *reader = vtkExtMisc::SafeReadPolyData(fname);
		if (!reader) return false;
		InitialiseSurface(reader->GetOutput());
		reader->Delete();
		ReadTexture(fname);
		return true;
	}
	else if (ext == "stl")
	{
		vtkSTLReader *reader = vtkExtMisc::SafeReadPolyDataSTL(fname);
		if (!reader) return false;
		InitialiseSurface(reader->GetOutput());
		reader->Delete();
		ReadTexture(fname);
		return true;
	}
	else if (ext == "obj")
	{
		vtkOBJReader *reader = vtkExtMisc::SafeReadPolyDataOBJ(fname);
		if (!reader) return false;
		InitialiseSurface(reader->GetOutput());
		reader->Delete();
		ReadTexture(fname);
		return true;
	}
	else if (ext == "wrl")
	{
		vtkVRMLImporter *vrmlin = vtkVRMLImporter::New();
		vrmlin->SetFileName(fname.c_str());
		vrmlin->Update();

		vtkPolyData *pd = (vtkPolyData*)(vrmlin->GetRenderer()->GetActors()->GetLastActor()->GetMapper()->GetInput());
		pd->GetPointData()->SetScalars(NULL); 
		InitialiseSurface(pd);
		vrmlin->Delete();
		ReadTexture(fname);
		return true;
	}
	else  // We do not use the multireader for stl, vtk, obj because of the texture
	{
		vtkPolyData *pd = vtkExtMisc::MultiReadSurface(fname);
		if (pd == NULL)
		{
			return false;
		}
		InitialiseSurface(pd);
		pd->Delete();
		return true;
	}
	return false;
}

bool CSurfaceProperties::file_exists( char const* fn )
{
	struct stat fs;
	return stat(fn, &fs) == 0;
}
