#ifndef _Landmark_h_
#define _Landmark_h_

#include <string>

//! Landmark class
/** Contains the same attributes as landmarks set using MeshLabs point picker 
 <point x="-56.7457" y="77.0165" z="-28.3349" active="1" name="1 Right brow start"/> */
class CLandmark
{
	public:
		//! Default constructor
		CLandmark();

		//! Clear pos and set to not active
		void Clear();

		//! The id number of the landmark (normally in which order it is placed). Could be seen as redundant info
		int id;

		//! x, y, z
		double pos[3];

		//! active (normally bool, but with int there are more options 
		int active;

		//! Name
		std::string name;
};

#endif
