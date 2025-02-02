/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkConstrainedICPTransform.h,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkConstrainedICPTransform - Implementation of the ICP algorithm.
// .SECTION Description
// Match two surfaces using the iterative closest point (ICP) algorithm.
// The core of the algorithm is to match each vertex in one surface with 
// the closest surface point on the other, then apply the transformation
// that modify one surface to best match the other (in a least square sense).
// This has to be iterated to get proper convergence of the surfaces.
// .SECTION Note
// Use vtkTransformPolyDataFilter to apply the resulting ICP transform to 
// your data. You might also set it to your actor's user transform.
// .SECTION Note
// This class makes use of vtkLandmarkTransform internally to compute the
// best fit. Use the GetLandmarkTransform member to get a pointer to that
// transform and set its parameters. You might, for example, constrain the
// number of degrees of freedom of the solution (i.e. rigid body, similarity,
// etc.) by checking the vtkLandmarkTransform documentation for its SetMode
// member.
// .SECTION see also
// vtkLandmarkTransform


#ifndef __vtkConstrainedICPTransform_h
#define __vtkConstrainedICPTransform_h

#include "vtkLinearTransform.h"
#include <vector>

#define VTK_ICP_MODE_RMS 0
#define VTK_ICP_MODE_AV 1

class vtkCellLocator;
class vtkLandmarkTransform;
class vtkDataSet;
class vtkPolyData;

class vtkConstrainedICPTransform : public vtkLinearTransform
{
public:
	static vtkConstrainedICPTransform *New();
	vtkTypeMacro(vtkConstrainedICPTransform,vtkLinearTransform);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	// Specify the source and target data sets.
	void SetSource(vtkDataSet *source);
	void SetTarget(vtkDataSet *target);
	vtkGetObjectMacro(Source, vtkDataSet);
	vtkGetObjectMacro(Target, vtkDataSet);

	//// Description:
	//// Set/Get a spatial locator for speeding up the search process. 
	//// An instance of vtkCellLocator is used by default.
	//void SetLocator(vtkCellLocator *locator);
	//vtkGetObjectMacro(Locator,vtkCellLocator);

	// Description: 
	// Set/Get the  maximum number of iterations
	vtkSetMacro(MaximumNumberOfIterations, int);
	vtkGetMacro(MaximumNumberOfIterations, int);

	// Description: 
	// Get the number of iterations since the last update
	vtkGetMacro(NumberOfIterations, int);

	// Description: 
	// Force the algorithm to check the mean distance between two iteration.
	vtkSetMacro(CheckMeanDistance, int);
	vtkGetMacro(CheckMeanDistance, int);
	vtkBooleanMacro(CheckMeanDistance, int);

	// Description:
	// Specify the mean distance mode. This mode expresses how the mean 
	// distance is computed. The RMS mode is the square root of the average
	// of the sum of squares of the closest point distances. The Absolute
	// Value mode is the mean of the sum of absolute values of the closest
	// point distances.
	vtkSetClampMacro(MeanDistanceMode,int,
		VTK_ICP_MODE_RMS,VTK_ICP_MODE_AV);
	vtkGetMacro(MeanDistanceMode,int);
	void SetMeanDistanceModeToRMS()
	{this->SetMeanDistanceMode(VTK_ICP_MODE_RMS);}
	void SetMeanDistanceModeToAbsoluteValue()
	{this->SetMeanDistanceMode(VTK_ICP_MODE_AV);}
	const char *GetMeanDistanceModeAsString();

	// Description: 
	// Set/Get the maximum mean distance between two iteration. If the mean
	// distance is lower than this, the convergence stops.
	vtkSetMacro(MaximumMeanDistance, double);
	vtkGetMacro(MaximumMeanDistance, double);

	// Description: 
	// Get the mean distance between the last two iterations.
	vtkGetMacro(MeanDistance, double);

	vtkSetMacro(DebugOutput, bool);


	// Description: 
	// Set/Get the maximum number of landmarks sampled in your dataset.
	// If your dataset is dense, then you will typically not need all the 
	// points to compute the ICP transform. 
	vtkSetMacro(MaximumNumberOfLandmarks, int);
	vtkGetMacro(MaximumNumberOfLandmarks, int);

	// 0 - cellLocator, 1 - pointlocator
	vtkSetMacro(LocatorType, int);
	vtkGetMacro(LocatorType, int);

	// 
	vtkSetMacro(DistanceThreshold, double);
	vtkGetMacro(DistanceThreshold, double);

	vtkSetMacro(CheckDistance, int);
	vtkGetMacro(CheckDistance, int);
	vtkBooleanMacro(CheckDistance, int);


	// Description: 
	// Starts the process by translating source centroid to target centroid.
	vtkSetMacro(StartByMatchingCentroids, int);
	vtkGetMacro(StartByMatchingCentroids, int);
	vtkBooleanMacro(StartByMatchingCentroids, int);

	// Description: 
	// Get the internal landmark transform. Use it to constrain the number of
	// degrees of freedom of the solution (i.e. rigid body, similarity, etc.).
	vtkGetObjectMacro(LandmarkTransform,vtkLandmarkTransform);

	vtkGetObjectMacro(Quarternion, double);
	vtkGetObjectMacro(Translation, double);

	vtkGetObjectMacro(CheckMatrix, vtkMatrix4x4);

	vtkGetMacro(MeanDistanceAfterThreshold, double);

	// Description:
	// Invert the transformation.  This is done by switching the
	// source and target.
	void Inverse();

	// Description:
	// Make another transform of the same type.
	vtkAbstractTransform *MakeTransform();

protected:

	// Description:
	// Release source and target
	void ReleaseSource(void);
	void ReleaseTarget(void);

	//// Description:
	//// Release locator
	//void ReleaseLocator(void);

	//// Description:
	//// Create default locator. Used to create one when none is specified.
	//void CreateDefaultLocator(void);

	// Description:
	// Get the MTime of this object also considering the locator.
	vtkMTimeType GetMTime();

	vtkConstrainedICPTransform();
	~vtkConstrainedICPTransform();

	void InternalUpdate();

	// Description:
	// This method does no type checking, use DeepCopy instead.
	void InternalDeepCopy(vtkAbstractTransform *transform);

	void CreateFeatureEdges();

	vtkDataSet* Source;
	vtkDataSet* Target;
//	vtkCellLocator *Locator;
	int MaximumNumberOfIterations;
	int CheckMeanDistance;
	int MeanDistanceMode;
	double MaximumMeanDistance;
	int MaximumNumberOfLandmarks;
	int StartByMatchingCentroids;

	// 0 - cellLocator, 1 - pointlocator
	int LocatorType;

	int NumberOfIterations;
	double MeanDistance;
	vtkLandmarkTransform *LandmarkTransform;

	vtkPolyData *featureEdges;

	vtkCellLocator *featureLocator;

	double Translation[3];
	double Quarternion[4];
	// Keeps track of borderpoints
	std::vector<int> BorderPointsID;

	// Test matrix for debug purposes
	vtkMatrix4x4 *CheckMatrix;

	bool DebugOutput;

	// Tolerance for finding edges
	double EdgeFinderTolerance;

	int CheckDistance;

	double DistanceThreshold;

	double MeanDistanceAfterThreshold;

private:
	vtkConstrainedICPTransform(const vtkConstrainedICPTransform&);  // Not implemented.
	void operator=(const vtkConstrainedICPTransform&);  // Not implemented.
};

#endif
