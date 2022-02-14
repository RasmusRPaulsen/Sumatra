#include "MultiTimer.h"

#undef min
#undef max

#include <iostream>
#include <algorithm>

CMultiTimer::CMultiTimer()
{
	InitializeCriticalSection(&m_cs);
	m_MaxTimes = 100;
	int RingBufferSize = 30;
	m_StartTimes.resize(m_MaxTimes, 0);
	m_EndTimes.resize(m_MaxTimes, 0);
	m_EventNames.resize(m_MaxTimes);
	m_ElapsedTimes.resize(m_MaxTimes);
	m_ETRingBuffer.resize(m_MaxTimes);
	m_ETRingBufferHead.resize(m_MaxTimes, 0);
	m_EventImportance.resize(m_MaxTimes, 1);
	m_NumMeas.resize(m_MaxTimes, 0);

	for (unsigned int i = 0; i < m_ETRingBuffer.size(); i++)
	{
		m_ETRingBuffer[i].resize(RingBufferSize);
	}

	m_idx = 0;
	m_lastID = -1;

	::QueryPerformanceFrequency( (LARGE_INTEGER*)&m_nFrequency );
}


CMultiTimer &CMultiTimer::MultiTimer()
{
	static CMultiTimer theMultiTimer;
	return theMultiTimer;
}


CMultiTimer::~CMultiTimer()
{
	DeleteCriticalSection(&m_cs);
}

void CMultiTimer::Start(const std::string& eventName, int eventImportance)
{
	EnterCriticalSection(&m_cs);
	int id = -1;
	for (int i = 0; i < m_idx && (id == -1); i++)
	{
		if (eventName == m_EventNames[i])
			id = i;
	}
	if (id == -1)
	{
		if (m_idx < m_MaxTimes-1)
		{
			id = m_idx++;
		}
	}
	if (id != -1)
	{
		double t;

		__int64 timeCurr;
		::QueryPerformanceCounter( (LARGE_INTEGER*)&timeCurr );

		t = (double)timeCurr / ((double)m_nFrequency);

		m_StartTimes[id] = t;
		m_EventNames[id] = eventName;
		m_EventImportance[id] = eventImportance;
		m_lastID = id;
	}
	LeaveCriticalSection(&m_cs);
}

void CMultiTimer::End(const std::string& eventName)
{
	EnterCriticalSection(&m_cs);
	int id = -1;
	if (m_lastID >= 0)
	{
		if (eventName == m_EventNames[m_lastID])
			id = m_lastID;
	}
	if (id == -1)
	{
		for (int i = 0; i < m_idx && (id == -1); i++)
		{
			if (eventName == m_EventNames[i])
				id = i;
		}
	}

	if (id != -1)
	{
		double t;

		__int64 timeCurr;
		::QueryPerformanceCounter( (LARGE_INTEGER*)&timeCurr );

		t = (double)timeCurr / ((double)m_nFrequency);

		m_EndTimes[id] = t;
		m_ElapsedTimes[id] = m_EndTimes[id] - m_StartTimes[id];

		m_ETRingBuffer[id][m_ETRingBufferHead[id]] = m_ElapsedTimes[id];
		if (++m_ETRingBufferHead[id] >= (int)m_ETRingBuffer[id].size())
		{
			m_ETRingBufferHead[id] = 0;
		}
		m_NumMeas[id]++;
	}
	LeaveCriticalSection(&m_cs);
}

void CMultiTimer::PrintAllTimes(std::ostream& ost)
{
	EnterCriticalSection(&m_cs);
	if (m_idx == 0)
	{
		std::cerr << "No recorded timings" << std::endl;
	}
	double minTime, maxTime;

	for (int i = 0; i < m_idx; i++)
	{
		int Importance = 0;
		GetEventStats(i, minTime, maxTime, Importance);
		ost << m_EventNames[i] << " " 
			<< m_ElapsedTimes[i] * 1000 
			<< " (" << minTime * 1000 << ", " << maxTime * 1000 << ") "
			<< " ms" << std::endl;
	}
	LeaveCriticalSection(&m_cs);
}


int CMultiTimer::GetNumberOfTimings() 
{
	EnterCriticalSection(&m_cs);
	int id = m_idx;
	LeaveCriticalSection(&m_cs);

	return id;
}

std::string CMultiTimer::GetEventName(int idx) 
{
	EnterCriticalSection(&m_cs);
	std::string en = "";
	if (idx >= 0 && idx < m_idx)
		en = m_EventNames[idx];
	LeaveCriticalSection(&m_cs);
	
	return en;
}

double CMultiTimer::GetEventTime(int idx) 
{
	EnterCriticalSection(&m_cs);
	double t = 0;
	if (idx >= 0 && idx < m_idx)
		t = m_ElapsedTimes[idx];
	LeaveCriticalSection(&m_cs);
	
	return t;
}

int CMultiTimer::GetEventImportance(int idx) 
{
	EnterCriticalSection(&m_cs);
	int imp = 0;
	if (idx >= 0 && idx < m_idx)
		imp = m_EventImportance[idx];
	LeaveCriticalSection(&m_cs);
	
	return imp;
}


void CMultiTimer::GetEventStats(int idx, double &minTime, double &maxTime, int &importance)
{
	EnterCriticalSection(&m_cs);
	if (idx >= 0 && idx < m_idx)
	{
		minTime = m_ETRingBuffer[idx][0];
		maxTime = m_ETRingBuffer[idx][0];

		for (unsigned int i = 0; i < std::min((unsigned int)m_ETRingBuffer[idx].size(), m_NumMeas[idx]); i++)
		{
			minTime = std::min(minTime, m_ETRingBuffer[idx][i]);
			maxTime = std::max(maxTime, m_ETRingBuffer[idx][i]);
		}
		importance = m_EventImportance[idx];
	}
	LeaveCriticalSection(&m_cs);
}

void CMultiTimer::Clear()
{
	EnterCriticalSection(&m_cs);
	m_lastID = -1;
	m_idx = 0;
	LeaveCriticalSection(&m_cs);
}


int CMultiTimer::GetNamedEventID(const std::string& eventName)
{
	EnterCriticalSection(&m_cs);
	int id = -1;
	if (m_lastID >= 0)
	{
		if (eventName == m_EventNames[m_lastID])
			id = m_lastID;
	}
	if (id == -1)
	{
		for (int i = 0; i < m_idx && (id == -1); i++)
		{
			if (eventName == m_EventNames[i])
				id = i;
		}
	}
	LeaveCriticalSection(&m_cs);
	return id;
}

