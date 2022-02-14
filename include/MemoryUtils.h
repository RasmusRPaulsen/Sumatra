#ifndef _MemoryUtils_h_
#define _MemoryUtils_h_

#define NOMINMAX
#include <windows.h>

//! Memory utilities 
class CMemoryUtils
{
	public:

		static void PrintMemoryStatus();

		static DWORDLONG GetAvailableMemoryForDouble();

};

#endif
