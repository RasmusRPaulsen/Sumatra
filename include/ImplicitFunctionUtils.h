#ifndef _ImplicitFunctionUtils_h_
#define _ImplicitFunctionUtils_h_

class vtkImplicitFunction;
class vtkPolyData;

//! Utilities used with implicit functions
/**  */
class CImplicitFunctionUtils
{
	public:
		//! Default constructor
		CImplicitFunctionUtils();

		//! Destructor
		virtual ~CImplicitFunctionUtils();

		//! Find the closest point (cp) to the point (p) on the zero level of the function
		/** Returns false in case of errors. Sanityjump is a threshold for how far a point should be able to move. Could for example be the closest point neighbour distance,*/
		static bool FindClosestPoint(vtkImplicitFunction *func, double *p, double *cp, double accuracy, double SanityJump, double MaxInitialDist = 1e8, int iterations = 1000);

		//! Find the closest point (cp) to the point (p) on the zero level of the function
		static bool FindClosestPointByBisection(vtkImplicitFunction* func, double* p, double* cp, double epsilon);

		//! Find the closest point (cp) to the point (p) on the zero level of the function
		/** Returns false in case of errors.
			The search path is returned */
		static bool FindClosestPointIllustrated(vtkImplicitFunction *func, double *p, double *cp, vtkPolyData *path, double accuracy = 0.1);

	private:
};

#endif
