#ifndef _ColorManager_h_
#define _ColorManager_h_

#include <vector>

//! Basic class
/** This is a class used as template for new classes */
class CColorManager
{
	public:
		//! Default constructor
		/** It tries to read a file called ColorManager.ini in the program directory.
		    If not succesfull it will return one very boring color always (white) */
		CColorManager();

		//! Get next color
		/** \param range1 if this is true, each component will be scaled to 0-1*/
		void GetNextColor(double *col, bool range1 = true);

	private:
		
		//! Extremely simple color class
		class CColManRGB
		{
		public:
			//! Red component
			double R;

			//! Green component
			double G;

			//! Blue component
			double B;
		};

		//! Read color file
		bool ReadIniFile();

		//! Vector of colors
		std::vector<CColManRGB> m_colors;
		
		//! Current color id
		unsigned int m_curId;
};

#endif
