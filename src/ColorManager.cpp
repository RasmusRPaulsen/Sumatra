#include "ColorManager.h"

#include <string>
#include <fstream>
#include "GeneralUtils.h"

CColorManager::CColorManager()
{
	m_curId = 0;
	ReadIniFile();
}

void CColorManager::GetNextColor(double *col, bool range1)
{
	if (m_colors.size() == 0)
	{
		col[0] = 255;
		col[1] = 255;
		col[2] = 255;
	}
	else
	{
		if (m_curId >= m_colors.size())
			m_curId = 0;

		col[0] = m_colors[m_curId].R;
		col[1] = m_colors[m_curId].G;
		col[2] = m_colors[m_curId].B;
	
		m_curId++;
	}

	if (range1)
	{
		col[0] /= 255.0;
		col[1] /= 255.0;
		col[2] /= 255.0;
	}

}


bool CColorManager::ReadIniFile()
{
	std::string colorFile = CGeneralUtils::GetPathFromFilename(__argv[0]) + "ColorManager.ini";

	std::ifstream fist(colorFile.c_str());

	if (!fist)
		return false;
	bool stop = false;
	do 
	{
		CColManRGB tRGB;
		fist >> tRGB.R >> tRGB.G >> tRGB.B;
		if (fist.fail())
		{
			stop = true;
		}
		else
		{
			m_colors.push_back(tRGB);
		}

	} while(!stop);
		
	if (m_colors.size() == 0)
		return false;

	return true;
}
