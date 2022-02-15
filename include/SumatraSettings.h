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

		//! Get next color
		/** \param range1 if this is true, each component will be scaled to 0-1*/
		void GetNextColor(double* col, bool range1 = true);

		// unsigned byte RGB
		vtkColor3ub mBackgroundColor;

		int mColorbarType = 0;

		int mPointSize = 1;

	private:
		bool ParseJSON(const QJsonObject& json);

		int mCurrentColor = 0;

		std::vector<vtkColor3ub> mColors;
};

#endif
