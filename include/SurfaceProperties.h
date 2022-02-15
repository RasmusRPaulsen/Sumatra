#ifndef _SurfaceProperties_h_
#define _SurfaceProperties_h_

#include <string>
#include <vector>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSTLReader.h>
#include <vtkOBJReader.h>
#include <vtkTexture.h>

#include "SumatraSettings.h"

class vtkTexture;

//! Describes surface properties
class CSurfaceProperties
{
	public:
		//! Constructor with lookup table and settings
		CSurfaceProperties(vtkLookupTable *lut, CSumatraSettings *settings);

		//! Destructor
		~CSurfaceProperties();

		//! Read a surface from file
		bool ReadFromFile(const std::string& fname);

		//! Save a surface
		bool SaveToFile(const std::string& fname, bool ApplyUserTransform = false, bool Binary = false,
			bool WriteNormals = false, bool WriteScalars = false);

		//! Create and set default look
		void InitialiseSurface(vtkPolyData *pd);
		
		//! Set a texture image
		void SetTexture(vtkImageData *tex);

		//! Call this function if scalar range have been changed or scalars have been added
		void UpdateScalarProperties();

		//! Use this function before any function that require an undo
		void StoreState();

		//! Restore the stored state
		void RestoreState();

		//! Short name (filename without path and extenstion if read)
		std::string m_shortname;

		//! Filename f.ex.
		std::string m_fullname;

		//! The surface
		vtkPolyData *m_polyData;

		//! The mapper
		vtkPolyDataMapper *m_mapper;

		//! Texture (if any)
		vtkTexture *m_Texture;

		//! The actor
		vtkActor *m_actor;

		//! Pointer to lookup table
		vtkLookupTable *m_lookup;

		//! Set the scalar range
		void SetScalarRange(double *range);

		//! Returns false if no scalars
		bool GetScalarRange(double *range);

		//! Does it have scalars
		bool HaveScalars() const;

		//! Is undo available
		bool UndoAvailable() const;
		
	private:

		//! Try to read texture
		/** Will search for png and bmp files */
		bool ReadTexture(const std::string& fname);

		//! Simple file check
		bool file_exists( char const* fn );

		//! Constructor
		CSurfaceProperties();
	
		CSumatraSettings* mSettings;

		//! Pointer to a colormanager
		//CColorManager *m_colorManager;

		//! Poly data copy for undo operations
		vtkPolyData *m_UndoPolyData;

		//! Undo matrix
		vtkMatrix4x4 *m_UndoMatrix;

		//! Undo name
		std::string m_UndoName1;
		std::string m_UndoName2;

		//! Is undo available
		bool m_UndoAvailable;

		//! Default point size - should be made into a general structure with more options
		// int m_DefPointSize;
};

#endif
