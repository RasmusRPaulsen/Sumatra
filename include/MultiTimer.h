#ifndef _MultiTimer_h_
#define _MultiTimer_h_

#pragma warning(disable : 4786)

#include <iostream>
#include <string>
#include <vector>
#include <map>
#define NOMINMAX
#include <windows.h>
#undef _WINDOWS_ // VERY VERY Dirty. To avoid fatal error. I need windows to define CRITICAL_SECTION

//! Multi timer singleton class
/** Used for keeping track of timings across a project.

   Singleton implementation (only one instance allowed). User should access all methods through CMultiTimer::MultiTimer().
   Usage:
     CMultiTimer::MultiTimer.Start("Eventname");
     CMultiTimer::MultiTimer.End("Eventname");

 */
class CMultiTimer
{
	public:

		//! The main interface method
		static CMultiTimer& MultiTimer();

		//! Destructor
		virtual ~CMultiTimer();

		//! Records the start time of a named event
		void Start(const std::string& eventName, int eventImportance = 1);

		//! Records the end time of a named event
		void End(const std::string& eventName);

		//! Print all recorded event and their timings
		void PrintAllTimes(std::ostream& ost);

		//! Get number of recorded event timings
		int GetNumberOfTimings();

		//! Get the ID of a named event
		/** return -1 if event not found */
		int GetNamedEventID(const std::string& eventName);

		//! Return the name of event idx
		std::string GetEventName(int idx);

		//! Get the elapsed time between start and end for the event idx
		double GetEventTime(int idx);

		//! Get the importance of a given event
		int GetEventImportance(int idx);

		//! Get the timing stats for the event idx
		void GetEventStats(int idx, double &minTime, double &maxTime, int &importance);

		//! Clear all recorded timings
		void Clear();

	private:
		//! Default constructor
		CMultiTimer();

		//! Current max id
		int m_idx;

		//! Start times
		std::vector<double> m_StartTimes;

		//! End times
		std::vector<double> m_EndTimes;

		//! Elapsed times
		std::vector<double> m_ElapsedTimes;

		//! Number of measurements
		std::vector<unsigned int> m_NumMeas;

		//! Ring buffers of elapsed times
		/** Use to calculate statistics */
		std::vector<std::vector<double> > m_ETRingBuffer;

		//! Head of ring buffers
		std::vector<int> m_ETRingBufferHead;

		//! Event names
		std::vector<std::string> m_EventNames;

		//! Give the event an importance (1 - high, 2 - less high, etc)
		std::vector<int> m_EventImportance;

		//! Max number of recordings
		int m_MaxTimes;

		//! Cached last id
		int m_lastID;

		//! Critical section object to make the class re-entrant
		CRITICAL_SECTION m_cs;

		__int64 m_nFrequency;
};

#endif
