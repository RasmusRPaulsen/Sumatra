#ifndef _SumatraSettings_h_
#define _SumatraSettings_h_

//! Basic class
/** This is a class used as template for new classes */
class CSumatraSettings
{
	public:
		//! Default constructor
		CSumatraSettings();

		//! Destructor
		virtual ~CSumatraSettings();

		bool ReadSettings();

	private:
};

#endif
