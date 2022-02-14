#ifndef _LandmarkCollection_h_
#define _LandmarkCollection_h_

#include "LandMark.h"
#include <vector>
#include <set>

class vtkPolyData;

//! A collection of landmarks
/**  */
class CLandmarkCollection
{
	public:
		//! Default constructor
		CLandmarkCollection();

		//! Destructor
		virtual ~CLandmarkCollection();

		void InsernextLandmark(const CLandmark& lm);

		//! set landmark at position i in the landmark vector (a landmark at position 0 would typically have id 1)
		bool SetLandmark(int i, const CLandmark& lm);

		//! Get landmark at position i
		CLandmark GetLandmark(int i) const;

		//! Number of landmarks
		size_t GetNumberOfLandmarks() const;

		//! Clear landmarks
		void Clear();

		//! Read collection from Jens Fagertun anno file
		bool ReadFromAnnoFile(const std::string &name);
	
		//! Read collection from MeshLab XML file
		bool ReadFromXLMFile(const std::string &name);

		//! Return a new instance of a VTK polydata with the landmarks as vertices. Non active landmarks has scalar value 0 an the active value 1
		vtkPolyData *GetAsVTKPolydata() const;

		//! Return a new instance of a VTK polydata with the landmarks as vertices but only the vertices in the id list. Non active landmarks has scalar value 0 an the active value 1
		vtkPolyData *GetSelectedIDsAsVTKPolydata(const std::vector<int>& ids) const;

		//! Returns the ids where the current collection and the given landmarkcollection have overlapping points
		std::vector<int> FindOverlapIDs(const CLandmarkCollection& LM2);
		
		//! Returns the ids where the current collection and the given landmarkcollection have overlapping points. Restrict IDs to the ones in the IDSet
		std::vector<int> FindOverlapIDs(const CLandmarkCollection& LM2, std::set<int>& IDset);
		
		
	private:

		//! The landmarks
		std::vector<CLandmark> m_Landmarks;
};

#endif
