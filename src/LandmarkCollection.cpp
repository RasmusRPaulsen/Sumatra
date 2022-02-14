#include "LandmarkCollection.h"

#include <vtkXMLDataParser.h>
#include <vtkXMLDataElement.h>

#include <vtkpolydata.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <algorithm>
#include "GeneralUtils.h"


CLandmarkCollection::CLandmarkCollection()
{
	Clear();
}

CLandmarkCollection::~CLandmarkCollection()
{
}

void CLandmarkCollection::InsernextLandmark( const CLandmark& lm )
{
	m_Landmarks.push_back(lm);
}

bool CLandmarkCollection::SetLandmark( int i, const CLandmark& lm )
{
	if (i < 0 || i >= m_Landmarks.size())
		return false;

	m_Landmarks[i] = lm;
	return true;
}

CLandmark CLandmarkCollection::GetLandmark( int i ) const
{
	if (i < 0 || i >= m_Landmarks.size())
		return CLandmark(); // Return dummy landmark

	return m_Landmarks[i];

}

size_t CLandmarkCollection::GetNumberOfLandmarks() const
{
	return m_Landmarks.size();
}

void CLandmarkCollection::Clear()
{
	m_Landmarks.clear();
}

bool CLandmarkCollection::ReadFromAnnoFile(const std::string &name)
{
	std::ifstream fist(name.c_str());
	if (!fist)
	{
		std::cerr << "Could not open " << name << " for reading" << std::endl;
		return false;
	}

	int idx = 0;
	bool stop = false;
	do
	{
		std::string str;
		std::getline(fist, str);


		std::vector<double> nums;
		int nsp = CGeneralUtils::GetNumbersFromString(str, nums);
		
		if (!fist.fail() && !fist.eof() && nsp == 2)
		{
			CLandmark newlm;
			newlm.id = idx;
			newlm.pos[0] = nums[0];
			newlm.pos[1] = nums[1];
			newlm.pos[2] = nums[2];
			newlm.active = true;
			newlm.name = "";

			m_Landmarks.push_back(newlm);
			idx++;
		}
		else
			stop = true;
	} while (!stop);

	if (idx < 1)
	{
		return false;
	}

	return true;
}

bool CLandmarkCollection::ReadFromXLMFile( const std::string &name )
{
	vtkXMLDataParser *XMLdparser = vtkXMLDataParser::New();
	XMLdparser->SetFileName(name.c_str());
	if (!XMLdparser->Parse())
	{
		std::cerr << "problem parsing " << name << std::endl;
		XMLdparser->Delete();
		return false;
	}

	m_Landmarks.clear();
	vtkXMLDataElement *rootEl = XMLdparser->GetRootElement();

	int idx = 1;
	for (int i = 0; i < rootEl->GetNumberOfNestedElements(); i++)
	{
		vtkXMLDataElement *nestEl = rootEl->GetNestedElement(i);

		if (strcmp(nestEl->GetName(), "point") == 0)
		{
			double x = 0;
			double y = 0;
			double z = 0;
			int active = 0;
			std::string PointName;

			PointName = std::string(nestEl->GetAttribute("name"));

			if (!nestEl->GetScalarAttribute("x", x))
			{
				std::cerr << "problem parsing x coord" << std::endl;
				XMLdparser->Delete();
				return false;
			}
			if (!nestEl->GetScalarAttribute("y", y))
			{
				std::cerr << "problem parsing y coord" << std::endl;
				XMLdparser->Delete();
				return false;
			}
			if (!nestEl->GetScalarAttribute("z", z))
			{
				std::cerr << "problem parsing z coord" << std::endl;
				XMLdparser->Delete();
				return false;
			}
			if (!nestEl->GetScalarAttribute("active", active))
			{
				std::cerr << "problem parsing active attribute" << std::endl;
				XMLdparser->Delete();
				return false;
			}
			CLandmark newlm;
			newlm.id = idx;
			newlm.pos[0] = x;
			newlm.pos[1] = y;
			newlm.pos[2] = z;
			newlm.active = active;
			newlm.name = PointName;

			m_Landmarks.push_back(newlm);
			idx++;
		}
	}
	XMLdparser->Delete();

//	std::cout << "Read " << m_Landmarks.size() << " landmarks" << std::endl;
	return true;
}

vtkPolyData * CLandmarkCollection::GetAsVTKPolydata() const
{
	if (m_Landmarks.size() == 0)
		return NULL;

	vtkPolyData *lms = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *scalars = vtkDoubleArray::New();

	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < m_Landmarks.size(); i++)
	{
		vtkIdType id = pts->InsertNextPoint(m_Landmarks[i].pos);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);
		double val = (double)m_Landmarks[i].active;
		scalars->InsertNextTuple(&val);
	}

	lms->SetPoints(pts);
	pts->Delete();
	lms->SetVerts(verts);
	verts->Delete();
	lms->GetPointData()->SetScalars(scalars);
	scalars->Delete();

	return lms;
}

vtkPolyData * CLandmarkCollection::GetSelectedIDsAsVTKPolydata(const std::vector<int>& ids) const
{
	if (m_Landmarks.size() == 0 || ids.size() == 0)
		return NULL;

	vtkPolyData *lms = vtkPolyData::New();
	vtkPoints *pts = vtkPoints::New();
	vtkCellArray *verts = vtkCellArray::New();
	vtkDoubleArray *scalars = vtkDoubleArray::New();

	scalars->SetNumberOfComponents(1);

	for (int i = 0; i < ids.size(); i++)
	{
		int Pid = ids[i];
		vtkIdType id = pts->InsertNextPoint(m_Landmarks[Pid].pos);
		verts->InsertNextCell(1);
		verts->InsertCellPoint(id);
		double val = (double)m_Landmarks[Pid].active;
		scalars->InsertNextTuple(&val);
	}

	lms->SetPoints(pts);
	pts->Delete();
	lms->SetVerts(verts);
	verts->Delete();
	lms->GetPointData()->SetScalars(scalars);
	scalars->Delete();

	return lms;

}

std::vector<int> CLandmarkCollection::FindOverlapIDs( const CLandmarkCollection& LM2 )
{
	std::vector<int> overlap;

	size_t nlm = std::min(m_Landmarks.size(), LM2.GetNumberOfLandmarks());

	for (unsigned int i = 0; i < nlm; i++)
	{
		if (m_Landmarks[i].active && LM2.GetLandmark(i).active)
		{
			overlap.push_back(i);
		}
	}
	
	return overlap;
}

std::vector<int> CLandmarkCollection::FindOverlapIDs(const CLandmarkCollection& LM2, std::set<int>& IDset)
{
	std::vector<int> overlap;

	size_t nlm = std::min(m_Landmarks.size(), LM2.GetNumberOfLandmarks());

	for (unsigned int i = 0; i < nlm; i++)
	{
		if (m_Landmarks[i].active && LM2.GetLandmark(i).active)
		{
			if (IDset.find(i) != IDset.end())
			{
				overlap.push_back(i);
			}
		}
	}

	return overlap;
}
