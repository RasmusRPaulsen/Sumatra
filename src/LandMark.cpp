#include "Landmark.h"

CLandmark::CLandmark()
{
	Clear();
}

void CLandmark::Clear()
{
	id = -1;
	pos[0] = 0;
	pos[1] = 0;
	pos[2] = 0;
	active = 0;
	name = "";
}

