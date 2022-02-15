#ifndef _SumatraSettings_h_
#define _SumatraSettings_h_

#include <qjsonobject.h>
#include <vtkColor.h>
#include <vector>

//! Sumatra settings
/** Used for parsing settings file */
class CSumatraSettings
{
	public:
		//! Default constructor
		CSumatraSettings();

		//! Destructor
		virtual ~CSumatraSettings();

		bool ReadSettings();

	private:
		bool ParseJSON(const QJsonObject& json);

		// unsigned byte RGB
		vtkColor3ub mBackgroundColor;

		std::vector<vtkColor3ub> mColors;
};

#endif
