#ifndef _SumatraSettings_h_
#define _SumatraSettings_h_

#include <qjsonobject.h>
#include <vtkColor.h>
#include <vector>
#include <qstring.h>

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

		double mBackgroundColor[3];

		int mColorbarType = 0;

		int mPointSize = 1;

	private:
		bool ParseJSON(const QJsonObject& json);

		int mCurrentColor = 0;

		std::vector<vtkColor3ub> mColors;

		QString mTest;
};

#endif
