#include "GeneralUtils.h"

#include <sstream>
#include <algorithm>
#include <vtkMath.h>
#include <cmath>
#include <utility>
#include <vtkTransform.h>
#include <vtkPlane.h>
#include <cctype> // for tolower
#include <locale>


void CGeneralUtils::GetXYZfromOffset(int *dim, int offset, int &x, int &y, int &z)
{
	z = offset / (dim[0]  * dim[1]);
	int rest = offset - z * dim[0]  * dim[1];
	y = rest / dim[0];
	x = rest - y * dim[0];
}

double CGeneralUtils::SumSquaredDifference(const std::vector<double>& x, const std::vector<double>& y)
{
	const int n = x.size();
	if (y.size() != n || n == 0)
	{
		//! \todo Issue error
		return 0;
	}

	double sumSQ = 0;
	for (int i = 0; i < n; i++)
	{
		double dif = x[i] - y[i];
		sumSQ += dif * dif;
	}
	return sumSQ;
}

double CGeneralUtils::AverageAbsDifference(const std::vector<double>& x, const std::vector<double>& y)
{
	const int n = x.size();
	if (y.size() != n || n == 0)
	{
		//! \todo Issue error
		return 0;
	}

	double sumDif = 0;
	for (int i = 0; i < n; i++)
	{
		sumDif += fabs(x[i] - y[i]);
	}
	sumDif /= n;
	return sumDif;
}



double CGeneralUtils::Correlation(const std::vector<double>& x, const std::vector<double>& y)
{
	const int n = x.size();
	if (y.size() != n || n == 0)
	{
		//! \todo Issue error
		return 0;
	}
	// find means
	double ax = 0;
	double ay = 0;
	for (int i = 0; i < n; i++)
	{
		ax += x[i];
		ay += y[i];
	}
	
	ax /= n;
	ay /= n;

	double sxx = 0;
	double syy = 0;
	double sxy = 0;
	for (int i = 0; i < n; i++)
	{
		double xt = x[i] - ax;
		double yt = y[i] - ay;
		sxx += xt*xt;
		syy += yt*yt;
		sxy += xt*yt;
	}
	if (sxx * syy == 0)
	{
		//! \todo Throw warning
		return 0;
	}
	double r = sxy / sqrt(sxx*syy);

	return r;
}


double CGeneralUtils::SDev(std::vector<double>& x)
{
	const int n = x.size();
	if (n == 0)
	{
		return 0;
	}
	double sx = 0;
	double sxx = 0;
	for (int i = 0; i < n; i++)
	{
		sx += x[i];
		sxx += x[i] * x[i];
	}
	double mean = sx / n;
	double sdev = sqrt(sxx/n - mean * mean);
	return sdev;
}

void CGeneralUtils::NoiseRobustMeanAndSdev(std::vector<double>& x, double lowFrac, double hiFrac, double &mean, double& sdev)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}

	int lowP = (int)(n * lowFrac);
	int hiP  = (int)(n * hiFrac);

	if (lowP > hiP || lowP < 0 || hiP >= n)
	{
		std::cerr << "fractiles incorrect" << std::endl;
		return;
	}

	std::sort(x.begin(), x.end());

	double sx = 0;
	double sxx = 0;
	for (int i = lowP; i <= hiP; i++)
	{
		sx += x[i];
		sxx += x[i] * x[i];
	}
	int N = (hiP-lowP) + 1;

	mean = sx / N;
	sdev = sqrt(sxx/N - mean * mean);
}

void CGeneralUtils::MeanAndSdev(const std::vector<double>& x, double &mean, double& sdev)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	MeanAndVariance(x, mean, sdev);
	sdev = sqrt(sdev);
}

void CGeneralUtils::MeanAndSdev(const std::vector<unsigned char>& x, double &mean, double& sdev)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	MeanAndVariance(x, mean, sdev);
	sdev = sqrt(sdev);
}

void CGeneralUtils::MeanAndSdev( const std::vector<unsigned char>& x, double &mean, double& sdev, int minx, int maxx )
{
	MeanAndVariance(x, mean, sdev, minx, maxx);
	sdev = sqrt(sdev);
}

void CGeneralUtils::MeanAndVariance(const std::vector<double>& x, double &mean, double& var)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	double sx = 0;
	double sxx = 0;
	for (int i = 0; i < n; i++)
	{
		sx += x[i];
		sxx += x[i] * x[i];
	}
	mean = sx / n;
	var = sxx/n - mean * mean;
}

void CGeneralUtils::MeanAndVariance(const std::vector<unsigned char>& x, double &mean, double& var)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	double sx = 0;
	double sxx = 0;
	for (int i = 0; i < n; i++)
	{
		sx += x[i];
		sxx += x[i] * x[i];
	}
	mean = sx / n;
	var = sxx/n - mean * mean;
}

void CGeneralUtils::MeanAndVariance( const std::vector<unsigned char>& x, double &mean, double& var, int minx, int maxx )
{
	const int n = maxx-minx;
	if (n <= 0)
	{
		return;
	}
	double sx = 0;
	double sxx = 0;
	for (int i = minx; i < maxx; i++)
	{
		sx += x[i];
		sxx += x[i] * x[i];
	}
	mean = sx / n;
	var = sxx/n - mean * mean;
}

void CGeneralUtils::VectorMinMax(std::vector<double>& x, double &minval, double& maxval)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	minval = x[0];
	maxval = x[0];
	for (unsigned int i = 0; i < x.size(); i++)
	{
		minval = std::min(x[i], minval);
		maxval = std::max(x[i], maxval);
	}
}


void CGeneralUtils::Median(std::vector<double>& x, double fractile, double& median)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	int med = (int)(n * fractile);
	if (med < 0 || med >= n)
	{
		std::cerr << "fractile incorrect" << std::endl;
		return;
	}

	std::sort(x.begin(), x.end());
	median = x[med];
}

void CGeneralUtils::MinMax(const std::vector<double>& x, double &minval, double& maxval)
{
	const int n = x.size();
	if (n == 0)
	{
		return;
	}
	minval = x[0];
	maxval = x[0];
	for (int i = 0; i < n; i++)
	{
		minval = std::min(minval, x[i]);
		maxval = std::max(maxval, x[i]);
	}
}


void CGeneralUtils::MakeAnglePositive(double& a)
{
	while (a < 0)
	{
		a += Pi() * 2;
	}
	while (a > Pi() * 2)
	{
		a -= Pi() * 2;
	}
}

void CGeneralUtils::StandardiseSphericalCoordinates(double& Phi, double& Theta, double& Psi)
{
	MakeAnglePositive(Phi);
	MakeAnglePositive(Psi);
	MakeAnglePositive(Theta);

	// \todo Check if this way of standardising spherical coordinates is correct
	if (Theta > Pi())
	{
		Theta = 2 * Pi() - Theta;
		Phi += Pi();
		MakeAnglePositive(Phi);
	}
}

double CGeneralUtils::GaussRandomNumber()
{
	double rsq = 0;
	double v1 = 0;
	double v2 = 0;
	do
	{
		// Pick two uniform numbers in the square extending from -1 to +1 in each direction
		v1 = vtkMath::Random(-1, 1);;
		v2 = vtkMath::Random(-1, 1);;
		// See if they are in the unit circle, and if they are not, try again.
		rsq = v1*v1 + v2*v2;
	} 
	while (rsq >= 1.0 || rsq == 0);
	double fac = sqrt(-2.0 * log(rsq)/rsq);
	// Make the Box-Muller transform to get a normal deviate.

	return v1 * fac;
}

// taken from http://www.codeproject.com/string/stringsplit.asp?select=328279&df=100&forumid=2167&exp=0
int CGeneralUtils::SplitString(const std::string& input, const std::string& delimiter, std::vector<std::string>& results)
{
	int iPos = 0;
	int newPos = -1;
	int sizeS2 = delimiter.size();
	int isize = input.size();
	
	std::vector<int> positions;
	
	newPos = input.find (delimiter, 0);
	
	if( newPos < 0 ) { return 0; }
	
	int numFound = 0;
	
	while( newPos > iPos )
	{
		numFound++;
		positions.push_back(newPos);
		iPos = newPos;
		newPos = input.find (delimiter, iPos+sizeS2+1);
	}
	
	for( unsigned int i=0; i <= positions.size(); i++ )
	{
		std::string s;

		int offset = 0;
		if( i == 0 ) 
		{ 
			s = input.substr( i, positions[i] ); 
		}
		else
		{
			offset = positions[i-1] + sizeS2;
		}
		
		if( offset < isize )
		{
			if( i == positions.size() )
			{
				s = input.substr(offset);
			}
			else if( i > 0 )
			{
				s = input.substr( positions[i-1] + sizeS2, 
					positions[i] - positions[i-1] - sizeS2 );
			}
		}
		if( s.size() > 0 )
		{
			results.push_back(s);
		}
	}
	return numFound;
}

int CGeneralUtils::SplitStringIntoNumbers( const std::string& input, const std::string& delimiter, std::vector<double>& results )
{
	std::vector<std::string> nums;
	int nsp = CGeneralUtils::SplitString(input, delimiter, nums);

	for (unsigned int i = 0; i < nums.size(); i++)
	{
		double n = atof(nums[i].c_str());
		results.push_back(n);
	}
	return nsp;
}


double CGeneralUtils::ClosestPointToPlane(vtkPolyData *pd, double p0[3], double n[3], double *cp, int &pid)
{
	double mindist = -1;
	for (int i = 0; i < pd->GetNumberOfPoints(); i++)
	{
		double p[3];
		pd->GetPoint(i, p);
		
		double dist = vtkPlane::DistanceToPlane(p, n, p0);

		if (mindist < 0 || dist < mindist)
		{
			mindist = dist;
			cp[0] = p[0]; cp[1] = p[1];cp[2] = p[2];
			pid = i;
		}
	}

	return mindist;
}

double CGeneralUtils::RMS(const std::vector<double>& x)
{
	int i;
	const int N = x.size();
	if (N < 1)
	{
		return 0;
	}
	double sumsq = 0;
	for (i = 0; i < N; i++)
	{
		sumsq += x[i] * x[i];
	}
	sumsq /= N;
	return sqrt(sumsq);
}


std::string CGeneralUtils::StripPathAndExtensionFromFilename(const std::string& fname)
{
	std::string nname = fname;

	// strip path
	int ix1 = fname.rfind("\\");
	int ix2 = fname.rfind("/");

	// Compare which of them is biggest (if both symbols are used in the path)
	int ix = std::max(ix1, ix2);

	if (ix != std::string::npos)
	{
		nname.assign(fname, ix+1, fname.size());
	}

	// strip extension
	ix = nname.rfind(".");
	if (ix != std::string::npos)
	{
		nname.assign(nname, 0, ix);
	}
	return nname;
}

std::string CGeneralUtils::StripPathFromFilename(const std::string& fname)
{
	std::string nname = fname;

	// strip path
	int ix = fname.rfind("\\");
	if (ix != std::string::npos)
	{
		nname.assign(fname, ix+1, fname.size());
	}
	return nname;
}

std::string CGeneralUtils::GetPathFromFilename(const std::string& fname)
{
	std::string nname = "";

	// strip path
	int ix = fname.rfind("\\");
	if (ix != std::string::npos)
	{
		nname.assign(fname, 0, ix+1);
	}
	ix = fname.rfind("/");
	if (ix != std::string::npos)
	{
		nname.assign(fname, 0, ix+1);
	}
	return nname;
}

std::string CGeneralUtils::StripExtensionFromFilename(const std::string& fname)
{
	std::string nname = fname;

	// strip extension
	int ix = fname.rfind(".");
	if (ix != std::string::npos)
	{
		nname.assign(fname, 0, ix);
	}
	return nname;
}

std::string CGeneralUtils::GetExtensionFromFilename(const std::string& fname)
{
	std::string nname = fname;

	int ix = fname.rfind(".");
	if (ix == std::string::npos)
		return "";

	nname.assign(fname, ix+1, fname.length());

	return ToLower(nname);
}


std::string CGeneralUtils::StripWhiteSpace(const std::string& str )
{
	std::string res;
	int len = str.length();
	int i = 0;
	while( i < len && isspace( (unsigned char)str[i] ) )
		i++;
  
	if ( i == len )
	{
		res = "";
		return res;
	}
  
	int j = len - 1;
	while( j >= 0 && isspace((unsigned char)str[j] ) )
		j--;
  
	res.assign( str.c_str() + i, j - i + 1 );
	return res;
}

std::string CGeneralUtils::StripDoubleGnyf(const std::string& str)
{
	std::string res;
	int len = str.length();
	int i = 0;
	while( i < len && str[i] == '\"')
		i++;
  
	if ( i == len )
	{
		res = "";
		return res;
	}
  
	int j = len - 1;
	while( j >= 0 && str[j] == '\"')
		j--;
  
	res.assign( str.c_str() + i, j - i + 1 );
	return res;
}

void CGeneralUtils::StripNewline( std::string& str, const std::string delimiter )
{
	FindReplace( str, "\n", delimiter );
}

std::string CGeneralUtils::StripNewlineInline( const std::string str, const std::string delimiter )
{
	std::string str2 = str;
	StripNewline( str2, delimiter );
	return str2;
}

void CGeneralUtils::FindReplace( std::string& str, const std::string find, const std::string delimiter )
{
	int pos = 0;
	pos = str.find_first_of( find, pos );
	while (pos>0)
	{
		str.replace( pos, 1, delimiter );
		pos = str.find_first_of( find, pos );
	}
}

std::string CGeneralUtils::ToLower(const std::string& str)
{
	std::string res = str;

	std::transform (res.begin(),res.end(), res.begin(), tolower);

	return res;
}

int CGeneralUtils::UniqueNumbers(std::vector<int> &v)
{
	int unique = 0;
	const int N = v.size();
	std::vector<bool> used(N, false);

	for (int i = 0; i < N; i++)
	{
		if (!used[i])
		{
			int val = v[i];
			unique++;

			for (int j = i; j < N; j++)
			{
				if (v[j] == val)
				{
					used[j] = true;
				}
			}
		}
	}
	return unique;
}

void CGeneralUtils::WindowedMedianFilter( const std::vector<double>& x, std::vector<double>& out, int size )
{
	out.resize(x.size());

	int sz2 = size/2;

	for (int i = 0; i < (int)x.size(); i++)
	{
		std::vector<double> med;

		for (int j = -sz2; j <= sz2; j++)
		{
			if (i+j >= 0 && i+j < (int)x.size())
			{
				med.push_back(x[i+j]);
			}
		}
		double m;
		Median(med, 0.50, m);
		
		out[i] = m;
	}
}

void CGeneralUtils::WindowedAverageFilter( const std::vector<double>& x, std::vector<double>& out, int size )
{
	out.resize(x.size());

	int sz2 = size/2;

	for (int i = 0; i < (int)x.size(); i++)
	{
		std::vector<double> avg;

		for (int j = -sz2; j <= sz2; j++)
		{
			if (i+j >= 0 && i+j < (int)x.size())
			{
				avg.push_back(x[i+j]);
			}
		}
		double a;
		double s;
		MeanAndSdev(avg, a, s);

		out[i] = a;
	}
}

int CGeneralUtils::GetNumbersFromString( const std::string& input, std::vector<double>& results )
{
	std::istringstream ss(input);   
	double number = 0;   
	int num = 0;
	bool stop = false;
	do
	{
		ss >> number;     

		if (ss.fail( ))    
		{
			stop = true;
		}
		else
		{
			results.push_back(number);
		}
	}
	while (!stop);
	return results.size() - 1;
}


std::vector<int> CGeneralUtils::bounds(int parts, int mem) 
{
	std::vector<int>bnd;
	int delta = mem / parts;
	int reminder = mem % parts;
	int N1 = 0, N2 = 0;
	bnd.push_back(N1);
	for (int i = 0; i < parts; ++i) 
	{
		N2 = N1 + delta;
		if (i == parts - 1)
			N2 += reminder;
		bnd.push_back(N2);
		N1 = N2;
	}
	return bnd;
}
