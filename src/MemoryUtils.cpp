#include "MemoryUtils.h"
#define NOMINMAX
#include <windows.h>
#include <stdio.h>
#include <tchar.h>

void CMemoryUtils::PrintMemoryStatus()
{
	const int DIV = 1024;
	const int WIDTH = 7;
	MEMORYSTATUSEX statex;

	statex.dwLength = sizeof(statex);

	GlobalMemoryStatusEx(&statex);

	_tprintf(TEXT("There is  %*ld percent of memory in use.\n"),
		WIDTH, statex.dwMemoryLoad);
	_tprintf(TEXT("There are %*I64d total KB of physical memory.\n"),
		WIDTH, statex.ullTotalPhys / DIV);
	_tprintf(TEXT("There are %*I64d free  KB of physical memory.\n"),
		WIDTH, statex.ullAvailPhys / DIV);
	_tprintf(TEXT("There are %*I64d total KB of paging file.\n"),
		WIDTH, statex.ullTotalPageFile / DIV);
	_tprintf(TEXT("There are %*I64d free  KB of paging file.\n"),
		WIDTH, statex.ullAvailPageFile / DIV);
	_tprintf(TEXT("There are %*I64d total KB of virtual memory.\n"),
		WIDTH, statex.ullTotalVirtual / DIV);
	_tprintf(TEXT("There are %*I64d free  KB of virtual memory.\n"),
		WIDTH, statex.ullAvailVirtual / DIV);

	// Show the amount of extended memory available.

	_tprintf(TEXT("There are %*I64d free  KB of extended memory.\n"),
		WIDTH, statex.ullAvailExtendedVirtual / DIV);
}

DWORDLONG CMemoryUtils::GetAvailableMemoryForDouble()
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);

	DWORDLONG doubleSize = sizeof(double);
	DWORDLONG doublefree = statex.ullAvailPhys / doubleSize;

	return doublefree;
}
