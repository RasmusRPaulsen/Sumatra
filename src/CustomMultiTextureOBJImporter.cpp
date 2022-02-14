#include "CustomMultiTextureOBJImporter.h"

#include <GeneralUtils.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkExtMisc.h>
#include <vtkPointData.h>
#include <vtkPNGReader.h>
#include <vtkJPEGReader.h>
#include <vtkTIFFReader.h>
#include <vtkBMPReader.h>
#include <vtkTexture.h>
#include <vtkImageData.h>


CCustomMultiTextureOBJImporter::CCustomMultiTextureOBJImporter()
{
}

CCustomMultiTextureOBJImporter::~CCustomMultiTextureOBJImporter()
{
	for (unsigned int i = 0; i < Meshes.size(); i++)
	{
		if (Meshes[i])
			Meshes[i]->Delete();
	}

	for (unsigned int i = 0; i < Textures.size(); i++)
	{
		if (Textures[i])
			Textures[i]->Delete();
	}
}

bool CheckForToken(const std::string &str, const std::string &tok, std::string &rest)
{
	std::string tt = CGeneralUtils::StripWhiteSpace(str);

	if (tt.size() < tok.size())
		return false;

	std::string tt2 = CGeneralUtils::ToLower(tt);
	size_t pos = tt2.find(tok);

	if (pos == 0) // Token has to be first word 
	{
		rest.assign(tt, tok.size(), tt.size());
		rest = CGeneralUtils::StripWhiteSpace(rest);
		return true;
	}
	return false;
}

bool CCustomMultiTextureOBJImporter::ReadFile( const std::string& fname )
{
// 	std::string MTLname = CGeneralUtils::StripExtensionFromFilename(fname) + ".mtl";
// 	if (!ReadMTLFile(MTLname))
// 		return false;
// 
// 	std::cout << "Read " << Materials.size() << " different materials";

	std::ifstream fist(fname.c_str());
	if (!fist)
	{
		LastError << "Could not read " << fname;
		return false;
	}

	// The complete set of points
	vtkPoints *pts = vtkPoints::New();

	// The complete set of tcoords
	vtkFloatArray *tcoords = vtkFloatArray::New();
	tcoords->SetNumberOfComponents(2);

	// The points for the current part
	vtkPoints *newpts = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();
	vtkFloatArray *newtcoords = vtkFloatArray::New();
	newtcoords->SetNumberOfComponents(2);

	// Look up table for the points (position at the org id points to the new id)
	std::vector<int> PtsLookup;

	// Look up table for the texture coords (position at the org id points to the new id)
	std::vector<int> TCLookup;

	int MTLCount = 0;
	bool stop = false;
	do 
	{
		std::string tt;
		std::getline(fist, tt);

		if (!fist.eof() && !fist.fail())
		{
			tt=CGeneralUtils::StripWhiteSpace(tt);

			// Material lib
			std::string rest;
			if (CheckForToken(tt, "mtllib", rest))
			{
				if (rest.size() < 1)
				{
					LastError << "Error parsing OBJ files. mtllib not followed by filename. Found " << rest;
					return false;
				}
				std::string MTLname = CGeneralUtils::GetPathFromFilename(fname) + rest;
				if (!ReadMTLFile(MTLname))
					return false;
			}
		

			std::vector<std::string> st;
			int nst = CGeneralUtils::SplitString(tt, " ", st);

			//// Material lib
			//if (nst > 0 && st[0] == "mtllib")
			//{
			//	std::string rest;
			//	CheckForToken(tt, "mtllib",rest);

			//	if (nst != 1)
			//	{
			//		LastError << "Error parsing OBJ files. mtllib not followed by single name";
			//		return false;
			//	}
			//	std::string MTLname = CGeneralUtils::GetPathFromFilename(fname) + st[1];
			//	if (!ReadMTLFile(MTLname))
			//		return false;
			//}

			// Vertices
			if (nst == 3 && st[0] == "v")
			{
				double p[3];
				p[0] = 	atof(st[1].c_str());
				p[1] = 	atof(st[2].c_str());
				p[2] = 	atof(st[3].c_str());
				pts->InsertNextPoint(p);
			}

			// Texture coords
			if (nst == 2 && st[0] == "vt")
			{
				double tc[2];
				tc[0] = 	atof(st[1].c_str());
				tc[1] = 	atof(st[2].c_str());
				tcoords->InsertNextTuple(tc);
			}

			// Faces
			if (nst == 3 && st[0] == "f")
			{
				// if faces starts without a "newmtl" - again a hack
				if (PtsLookup.size() != pts->GetNumberOfPoints() || TCLookup.size() != tcoords->GetNumberOfTuples())
				{
					// Clear lookup table
					PtsLookup.resize(pts->GetNumberOfPoints(), -1);
					std::fill(PtsLookup.begin(), PtsLookup.end(), -1);   
					TCLookup.resize(tcoords->GetNumberOfTuples(), -1);
					std::fill(TCLookup.begin(), TCLookup.end(), -1);   

				}

				vtkIdType fi[3];

				for (int i = 0; i < 3; i++)
				{
					std::vector<std::string> st1;
					int nst1 = CGeneralUtils::SplitString(st[1+i], "//", st1);

					// Pointid//normalid
					if (nst1 == 1)
					{
						int idx = atoi(st1[0].c_str()) - 1; // subtract 1 due to 1-indexing in obj

						// lets check if this point id has already been used
						if (PtsLookup[idx] == -1)
						{
							 // add the point and put it in the lookup
							int newid = newpts->InsertNextPoint(pts->GetPoint(idx));
							PtsLookup[idx] = newid;
							
						}
						fi[i] = PtsLookup[idx];
					}
					else
					{
						st1.clear();
						nst1 = CGeneralUtils::SplitString(st[1+i], "/", st1);

						// Pointid/textureid
						if (nst1 == 1)
						{
							int idx = atoi(st1[0].c_str()) - 1; // subtract 1 due to 1-indexing in obj

							// lets check if this point id has already been used
							if (PtsLookup[idx] == -1)
							{
								// add the point and put it in the lookup
								int newid = newpts->InsertNextPoint(pts->GetPoint(idx));
								PtsLookup[idx] = newid;

							}
							fi[i] = PtsLookup[idx];

							// texture id
							idx = atoi(st1[1].c_str()) - 1; // subtract 1 due to 1-indexing in obj

							// lets check if this tc id has already been used
							if (TCLookup[idx] == -1)
							{
								// add the texture coords point and put it in the lookup
								int newid = newtcoords->InsertNextTuple(tcoords->GetTuple(idx));
								TCLookup[idx] = newid;
							}
						}
						// Pointid/textureid/normalid
						else if (nst1 == 2)
						{
							int idx = atoi(st1[0].c_str()) - 1; // subtract 1 due to 1-indexing in obj

							// lets check if this point id has already been used
							if (PtsLookup[idx] == -1)
							{
								// add the point and put it in the lookup
								int newid = newpts->InsertNextPoint(pts->GetPoint(idx));
								PtsLookup[idx] = newid;

							}
							fi[i] = PtsLookup[idx];

							// texture id
							idx = atoi(st1[1].c_str()) - 1; // subtract 1 due to 1-indexing in obj

							// lets check if this tc id has already been used
							if (TCLookup[idx] == -1)
							{
								// add the texture coords point and put it in the lookup
								int newid = newtcoords->InsertNextTuple(tcoords->GetTuple(idx));
								TCLookup[idx] = newid;
							}

							//if (st1[0] != st1[1])
							//	std::cerr << "pointid and textureid do not match " <<st1[0] << " - " << st1[1] << std::endl;
						}
						else // just the point
						{
							int idx = atoi(st[i+1].c_str()) - 1; // subtract 1 due to 1-indexing in obj

							// lets check if this point id has already been used
							if (PtsLookup[idx] == -1)
							{
								// add the point and put it in the lookup
								int newid = newpts->InsertNextPoint(pts->GetPoint(idx));
								PtsLookup[idx] = newid;

							}
							fi[i] = PtsLookup[idx];
						}
					}
				}

				polys->InsertNextCell(3,fi);
			}


			// usemtl
			if (nst == 1 && st[0] == "usemtl")
			{
				if (polys->GetNumberOfCells() > 0 && newpts->GetNumberOfPoints() > 0)
				{
					// First save old poly
					vtkPolyData *pd = vtkPolyData::New();
					pd->SetPoints(newpts);
					pd->SetPolys(polys);
					if (newtcoords->GetNumberOfTuples() > 0)
					{
						pd->GetPointData()->SetTCoords(newtcoords);
					}

					//std::ostringstream ost;
					//ost << "E:\\data\\IMM\\Vectra\\TestOut" << MTLCount << ".vtk";

					//vtkExtMisc::WritePDVTK(pd, ost.str());

					//std::cout << "Wrote " << ost.str() << std::endl << "with " << newpts->GetNumberOfPoints() << " points and " << polys->GetNumberOfCells() << " cells and " << newtcoords->GetNumberOfTuples() << " TC" << std::endl;
//					std::cout << "Created surface with " << newpts->GetNumberOfPoints() << " points and " << polys->GetNumberOfCells() << " cells and " << newtcoords->GetNumberOfTuples() << " TC" << std::endl;

					
					vtkPolyData *pd2 = vtkPolyData::New();
					pd2->DeepCopy(pd);
					Meshes.push_back(pd2);
					pd->Delete();

					MTLCount++;
					polys->Reset();
					polys->InitTraversal();
					newtcoords->Initialize();
					newpts->Initialize();
				}
				// Clear lookup table
				PtsLookup.resize(pts->GetNumberOfPoints(), -1);
				std::fill(PtsLookup.begin(), PtsLookup.end(), -1);   
				TCLookup.resize(tcoords->GetNumberOfTuples(), -1);
				std::fill(TCLookup.begin(), TCLookup.end(), -1);   
			}

			
			// Ignore all other tags
		}
		else
		{
			// Write last polys
			if (polys->GetNumberOfCells() > 0 && newpts->GetNumberOfPoints() > 0)
			{
				// First save old poly
				vtkPolyData *pd = vtkPolyData::New();
				pd->SetPoints(newpts);
				pd->SetPolys(polys);
				if (tcoords->GetNumberOfTuples() > 0)
				{
					pd->GetPointData()->SetTCoords(tcoords);
				}
				if (newtcoords->GetNumberOfTuples() > 0)
				{
					pd->GetPointData()->SetTCoords(newtcoords);
				}

				vtkPolyData *pd2 = vtkPolyData::New();
				pd2->DeepCopy(pd);
				Meshes.push_back(pd2);
				pd->Delete();

//				std::cout << "Created surface with " << newpts->GetNumberOfPoints() << " points and " << polys->GetNumberOfCells() << " cells and " << newtcoords->GetNumberOfTuples() << " TC" << std::endl;
				
				//std::ostringstream ost;
				//ost << "E:\\data\\IMM\\Vectra\\TestOut" << MTLCount << ".vtk";

				//vtkExtMisc::WritePDVTK(pd, ost.str());

				//std::cout << "Wrote " << ost.str() << std::endl << "with " << newpts->GetNumberOfPoints() << " points and " << polys->GetNumberOfCells() << " cells and " << newtcoords->GetNumberOfTuples() << " TC" << std::endl;

				//pd->Delete();
				MTLCount++;
				polys->Reset();
				polys->InitTraversal();
				// Clear lookup table
				PtsLookup.resize(pts->GetNumberOfPoints(), -1);
				std::fill(PtsLookup.begin(), PtsLookup.end(), -1);   
				TCLookup.resize(tcoords->GetNumberOfTuples(), -1);
				std::fill(TCLookup.begin(), TCLookup.end(), -1);   
			}
			stop = true;	
		}
	} while (!stop);

	pts->Delete();
	tcoords->Delete();
	polys->Delete();
	newtcoords->Delete();
	newpts->Delete();

	if (Meshes.size() == 0)
		return false;

	ReadTextures(fname);

	return true;
}

bool CheckIfImageFile(const std::string& fname)
{
	std::string ext = CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(fname));
	if (ext == "png" || ext == "jpg" || ext == "jpeg" || ext == "bmp" || ext == "tif" || ext == "tiff")
		return true;

	return false;
}

bool CCustomMultiTextureOBJImporter::ReadMTLFile( const std::string& fname )
{
	std::ifstream fist(fname.c_str());
	if (!fist)
	{
		LastError << "Could not read " << fname;
		return false;
	}

	// Simple check - we need more than 0
	int NImageFiles = 0;
	bool stop = false;
	do 
	{
		std::string tt;
		std::getline(fist, tt);
		if (!fist.eof() && !fist.fail())
		{
			tt=CGeneralUtils::StripWhiteSpace(tt);

			// New material
			std::string rest;
			if (CheckForToken(tt, "newmtl", rest))
			{
				if (rest.size() < 1)
				{
					LastError << "Error parsing MTL files. newmtl not followed by name. Found " << rest;
					return false;
				}
				Materials.push_back(MTLType(rest));
			}

			// Texture file name
			if (CheckForToken(tt, "map_kd", rest))
			{
				if (rest.size() < 1)
				{
					LastError << "Error parsing MTL files. map_Kd not followed by name. Found " << rest;
					return false;
				}
				if (Materials.size() == 0)
				{
					LastError << "newmtl needed before map_kd";
					return false;
				}
				if (!CheckIfImageFile(rest))
				{
					LastError << "This reader ONLY support map_kd followed by imagefile name. Found " << rest;
					return false;
				}

				NImageFiles++;
				Materials[Materials.size()-1].filename = rest;
			}

			//std::vector<std::string> st;
			//int nst = CGeneralUtils::SplitString(tt, " ", st);

			//// New material
			////if (nst > 0 && st[0] == "newmtl")
			////{
			////	if (nst != 1)
			////	{
			////		LastError << "Error parsing MTL files. newmtl not followed by single name";
			////		return false;
			////	}
			////	Materials.push_back(MTLType(st[1]));
			////}

			//// Texture file name
			//if (nst > 0 && st[0] == "map_Kd")
			//{
			//	if (nst != 1)
			//	{
			//		LastError << "Error parsing MTL files. map_Kd not followed by single name";
			//		return false;
			//	}
			//	if (Materials.size() == 0)
			//	{
			//		LastError << "newmtl needed before map_kd";
			//		return false;
			//	}
			//	if (!CheckIfImageFile(st[1]))
			//	{
			//		LastError << "This reader ONLY support map_kd followed by imagefile name";
			//		return false;
			//	}
			//	
			//	NImageFiles++;
			//	Materials[Materials.size()-1].filename = st[1];
			//}

			// Ignore all other tags
		}
		else
			stop = true;	
	} while (!stop);

	if (NImageFiles == 0)
	{
		LastError << "No image file names found with tag map_Kd";
		return false;
	}

	return true;
}


std::string CCustomMultiTextureOBJImporter::GetError() const
{
	return LastError.str();
}

vtkPolyData * CCustomMultiTextureOBJImporter::GetMesh( int i )
{
	if (i < 0 || i >= (int)Meshes.size())
		return NULL;
	return Meshes[i];
}

size_t CCustomMultiTextureOBJImporter::NumberOfMeshes() const
{
	return Meshes.size();
}


bool CCustomMultiTextureOBJImporter::ReadTextures( const std::string &fname )
{
	for (unsigned int i = 0; i < Materials.size(); i++)
	{
		vtkImageData *texture = NULL;
		std::string texname = Materials[i].filename;
		if (texname != "" && CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(Materials[i].filename)) == "png")
		{
			std::string TexFile = CGeneralUtils::GetPathFromFilename(fname) + texname;

			vtkPNGReader *textureImage = vtkPNGReader::New();
			textureImage->SetFileName(TexFile.c_str());
			textureImage->Update();
			
			if (textureImage)
			{
				texture = vtkImageData::New();
				texture->DeepCopy(textureImage->GetOutput());
			}
			textureImage->Delete();
		}
		else if (texname != "" && CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(Materials[i].filename)) == "jpg")
		{
			std::string TexFile = CGeneralUtils::GetPathFromFilename(fname) + texname;

			vtkJPEGReader *textureImage = vtkJPEGReader::New();
			textureImage->SetFileName(TexFile.c_str());
			textureImage->Update();

			if (textureImage)
			{
				texture = vtkImageData::New();
				texture->DeepCopy(textureImage->GetOutput());
			}
			textureImage->Delete();
		}
		else if (texname != "" && CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(Materials[i].filename)) == "tif" || CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(Materials[i].filename)) == "tiff")
		{
			std::string TexFile = CGeneralUtils::GetPathFromFilename(fname) + texname;

			vtkTIFFReader *textureImage = vtkTIFFReader::New();
			textureImage->SetFileName(TexFile.c_str());
			textureImage->Update();

			if (textureImage)
			{
				texture = vtkImageData::New();
				texture->DeepCopy(textureImage->GetOutput());
			}
			textureImage->Delete();
		}
		else if (texname != "" && CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(Materials[i].filename)) == "bmp" || CGeneralUtils::ToLower(CGeneralUtils::GetExtensionFromFilename(Materials[i].filename)) == "bmp")
		{
			std::string TexFile = CGeneralUtils::GetPathFromFilename(fname) + texname;

			vtkBMPReader *textureImage = vtkBMPReader::New();
			textureImage->SetFileName(TexFile.c_str());
			textureImage->Update();

			if (textureImage)
			{
				texture = vtkImageData::New();
				texture->DeepCopy(textureImage->GetOutput());
			}
			textureImage->Delete();
		}

		//if (texture)
		//{
		//	int* dims = texture->GetDimensions();
		//	std::cout << "Texture Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
		//	std::cout << "Data dim " << texture->GetDataDimension() << " number of com " << texture->GetNumberOfScalarComponents() << std::endl;

		//	int inc[3];
		//	texture->GetIncrements(inc);
		//	std::cout << "Increments " << inc[0] << ", " << inc[1] << ", " << inc[2] << std::endl;
		//	std::cout << "Data type " << texture->GetScalarTypeAsString() << std::endl;
		//}

		Textures.push_back(texture);
	}
	return true;
}

vtkImageData * CCustomMultiTextureOBJImporter::GetTexture( int i )
{
	if (i < 0 || i >= (int)Textures.size())
		return NULL;
	return Textures[i];

}

CCustomMultiTextureOBJImporter::MTLType::MTLType( const std::string &n )
{
	name = n;
}
