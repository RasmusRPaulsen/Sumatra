#ifndef _SumatraSettings_h_
#define _SumatraSettings_h_

#include <qjsonobject.h>

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

};

#endif
