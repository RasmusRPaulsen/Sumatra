#ifndef _GeneralUtils_h_
#define _GeneralUtils_h_

#include <vector>
#include <algorithm>
#include <functional>

#include <vtkPolyData.h>
#include <vtkMatrix4x4.h>

//! General utilities that only depends on (but minimally) 
class CGeneralUtils
{
	public:

		//! Useful constant
		static double DegreesToRadians() {return 0.017453292519943295;};

		//! Useful constant
		static double Pi() {return 3.1415926535897932384626;};

		//! Useful constant
		static double RadiansToDegrees() {return 57.29577951308232;};

		//! Electromagnetic constant
		static double VacuumPermeability() {return 4 * Pi() * 1e-7;};

		//! Calculate the average absolute difference between two series
		static double AverageAbsDifference(const std::vector<double>& x, const std::vector<double>& y);

		//! Calculate the correlation between two measurements series
		static double Correlation(const std::vector<double>& x, const std::vector<double>& y);

		//! Calculate the sum of squared differences between the two series
		static double SumSquaredDifference(const std::vector<double>& x, const std::vector<double>& y);

		//! Standard deviation
		static double SDev(std::vector<double>& x);

		//! Mean and standard deviation of a vector of doubles
		static void MeanAndSdev(const std::vector<double>& x, double &mean, double& sdev);
		
		//! Mean and standard deviation of a vector of bytes
		static void MeanAndSdev(const std::vector<unsigned char>& x, double &mean, double& sdev);

		//! Mean and standard deviation of a vector of bytes in a range
		static void MeanAndSdev(const std::vector<unsigned char>& x, double &mean, double& sdev, int minx, int maxx);

		//! Mean and variance
		static void MeanAndVariance(const std::vector<double>& x, double &mean, double& var);
		static void MeanAndVariance(const std::vector<unsigned char>& x, double &mean, double& var, int minx, int maxx);
		static void MeanAndVariance(const std::vector<unsigned char>& x, double &mean, double& var);
	
		//! Noise robust mean and standard deviation of a vector of doubles
		/** Sorts the vector and calculates the mean and sdev between
		    the low and high fractile.
		    Returns the sorted vector.*/
		static void NoiseRobustMeanAndSdev(std::vector<double>& x, double lowFrac, double hiFrac, double &mean, double& sdev);
		
		//! Median a vector of doubles
		/** Note: returns x as a sorted array*/
		static void Median(std::vector<double>& x, double fractile, double& median);

		//! Filter a vector using a local windowed median filter.
		/** Slow version */
		static void WindowedMedianFilter(const std::vector<double>& x, std::vector<double>& out, int size);

		//! Filter a vector using a local windowed average filter.
		/** Slow version */
		static void WindowedAverageFilter(const std::vector<double>& x, std::vector<double>& out, int size);

		//! Returns the Root Mean Square value
		static double RMS(const std::vector<double>& x);

		//! Min and max value of a vector of doubles
		static void MinMax(const std::vector<double>& x, double &minval, double& maxval);

		//! Put angle in 0 -> 2PI. Works in radians
		static void MakeAnglePositive(double& a);

		//! Standardise Spherical Coordinates
		/** After 
		    0 < Phi  < 2Pi
			0 < Theta < Pi
			0 < Psi < 2Pi. 
			Works in radians*/
		static void StandardiseSphericalCoordinates(double& Phi, double& Theta, double& Psi);

		//! Sort two vectors using so x will end up being in ascending order
		/** It is assumed that the vectors consist of (x,y) pairs */
		template <class T1, class T2>
		static void Sort2Vectors(std::vector<T1>& x, std::vector<T2>& y);		

		//! Sort two vectors using so x will end up being in descending order
		/** It is assumed that the vectors consist of (x,y) pairs */
		template <class T1, class T2>
		static void Sort2VectorsReverse(std::vector<T1>& x, std::vector<T2>& y);		

		//! Adapted from Numerical Recipes in C
		/** Returns random number with sdev = 1. */
		static double GaussRandomNumber();


		//! Finds min and max values in a vector of doubles
		static void VectorMinMax(std::vector<double>& x, double &minval, double& maxval);

		//! Splits a string into substrings based on a given delimter
		/** Can be used to parsing comma seperated files for example.
		    Returns the amount of delimeters found in the input string.*/
		static int SplitString(const std::string& input, const std::string& delimiter, std::vector<std::string>& results);
		
		//! Splits a string into substrings based on a given delimter
		/** Assumes it is made out of numbers */
		static int SplitStringIntoNumbers(const std::string& input, const std::string& delimiter, std::vector<double>& results);
		
		//! Given a string with numbers seperated by whitespace the numbers are returned
		/** Returns the amount of numbers-1 (==number of delimiters */
		static int GetNumbersFromString(const std::string& input, std::vector<double>& results);
				
		//! Return the point in a polydata that are closest to the plane.
		/** The plane is defined by n(x-p0) = 0. The normal n[3] must be magnitude=1*/
		static double ClosestPointToPlane(vtkPolyData *pd, double p0[3], double n[3], double *cp, int &pid);

		//! Strip path and extension from a file name
		static std::string StripPathAndExtensionFromFilename(const std::string& fname);

		//! Strip extension from a file name
		static std::string StripExtensionFromFilename(const std::string& fname);

		//! Strip path from a file name
		static std::string StripPathFromFilename(const std::string& fname);

		//! Get the path from a filename with ending "\"
		static std::string GetPathFromFilename(const std::string& fname);

		//! Get extension from file name (lowercased)
		static std::string GetExtensionFromFilename(const std::string& fname);

		//! Remove whitespaces from the start and end of a string
		static std::string StripWhiteSpace(const std::string& str);

		//! Remove " from the start and end of a string
		static std::string StripDoubleGnyf(const std::string& str);

		//! Replace newlines in the string by another delimiter
		static std::string StripNewlineInline( const std::string str, const std::string delimiter = " // " );
		static void StripNewline( std::string& str, const std::string delimiter = " // " );
		
		//! Find a substring and replace it by another string
		/** \todo This FindReplace seems broken */
		static void FindReplace( std::string& str, const std::string find, const std::string delimiter );

		//! Convert an entire string to lower-case
		static std::string ToLower(const std::string& str);

		//! Return the number of unique numbers in a vector
		/** Should be templated and made faster */
		static int UniqueNumbers(std::vector<int> &v);

		//! Split mem into parts
		/** e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10
		   if possible the function will split mem into equal chunks, if not 
		   the last chunk will be slightly larger */
		static std::vector<int> bounds(int parts, int mem);
		
		//!  Get coordinates (x, y, z) in a volume based on offset value
		static void GetXYZfromOffset(int *dim, int offset, int &x, int &y, int &z);
};

template <class T1, class T2>
void CGeneralUtils::Sort2Vectors(std::vector<T1>& x, std::vector<T2>& y)
{
	std::vector<std::pair<double,double> > t;

	for (unsigned int i = 0; i < x.size(); i++)
	{
		t.push_back(std::make_pair(x[i], y[i]));
	}
	std::sort(t.begin(), t.end());

	for (unsigned int i = 0; i < t.size(); i++)
	{
		x[i] = t[i].first;
		y[i] = t[i].second;
	}
}

template <class T1, class T2>
void CGeneralUtils::Sort2VectorsReverse(std::vector<T1>& x, std::vector<T2>& y)
{
	std::vector<std::pair<double,double> > t;

	for (unsigned int i = 0; i < x.size(); i++)
	{
		t.push_back(std::make_pair(x[i], y[i]));
	}
	std::sort(t.begin(), t.end(), std::greater<std::pair<double,double> > ( ));

	for (unsigned int i = 0; i < t.size(); i++)
	{
		x[i] = t[i].first;
		y[i] = t[i].second;
	}
}

#endif
