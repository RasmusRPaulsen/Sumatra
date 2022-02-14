#ifndef _CustomMultiTextureOBJImporter_h_
#define _CustomMultiTextureOBJImporter_h_

#include <string>
#include <vector>
#include <sstream>
class vtkPolyData;
class vtkImageData;

//! A custom importer of OBJ files that can also split a file into seperate polydatas
/** Currently customised to read Vectra M3 scans*/
class CCustomMultiTextureOBJImporter
{
	public:
		//! Default constructor
		CCustomMultiTextureOBJImporter();

		//! Destructor
		virtual ~CCustomMultiTextureOBJImporter();

		//! Give the OBJ file name. It will also search for the mtl file.
		bool ReadFile(const std::string& fname);

		//! get error
		std::string GetError() const;

		//! Get number of meshes
		size_t NumberOfMeshes() const;

		//! Return pointer to mesh
		vtkPolyData *GetMesh(int i);

		vtkImageData *GetTexture(int i);

	private:

		//! Very simple material type
		class MTLType
		{
		public:
			MTLType(const std::string &n);

			std::string name;
			std::string filename;
		};

		//! read material def file
		bool ReadMTLFile(const std::string& fname);

		//! Read textures
		bool ReadTextures(const std::string &fname);

		//! Vector of materials
		std::vector<MTLType> Materials;

		//! Keep track of error
		std::ostringstream LastError;

		//! The read meshes
		std::vector<vtkPolyData*> Meshes;

		//! The textures - not all meshes have a texture (will be NULL in the list then)
		std::vector<vtkImageData*> Textures;
};

#endif
