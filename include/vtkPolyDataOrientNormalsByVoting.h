// .NAME vtkPolyDataOrientNormalsByVoting - Make normals point the same direction using a graph based voting system
// .SECTION Description
// The point set is divided into sets of connected components (cliques) based on their inter-point distances
// and their mutual angle.
// It is possible to keep only the largest clique, cliques above a certain size, or all cliques
//
// It is also possible to supply a reference set of "normals" that will be used in the voting system.
// This could for example be a set of approximate normals given by a capture device.

#ifndef __vtkPolyDataOrientNormalsByVoting_h
#define __vtkPolyDataOrientNormalsByVoting_h

#include "vtkPolyDataAlgorithm.h"
#include <vector>

class vtkPolyData;

class  vtkPolyDataOrientNormalsByVoting : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPolyDataOrientNormalsByVoting,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	static vtkPolyDataOrientNormalsByVoting *New();

	// Description: 
	// The radius used when searching for neighbours
	vtkGetMacro(SearchRadius,double);
	vtkSetMacro(SearchRadius,double);

	// Description: 
	// The neighbour angle
	vtkGetMacro(NeighbourAngle,double);
	vtkSetMacro(NeighbourAngle,double);

	// Description: 
	// The minimum clique that will be kept
	// If equals zero all points are kept
	vtkGetMacro(MinCliqueSize,int);
	vtkSetMacro(MinCliqueSize,int);

	// Description:
	// If this is set to true, only the largest (most points) will be kept
	// no matter what MinCliqueSize is set to
	vtkSetMacro(LargestCliqueOnly,int);
	vtkGetMacro(LargestCliqueOnly,int);
	vtkBooleanMacro(LargestCliqueOnly,int);

	// Description:
	// This can be set and will be used to orient the normals
	// if NULL it is obviously not used
	vtkSetMacro(ReferencePolyData,vtkPolyData*);
	vtkGetMacro(ReferencePolyData,vtkPolyData*);

protected:
	vtkPolyDataOrientNormalsByVoting();
	~vtkPolyDataOrientNormalsByVoting() {};

	// Usual data generation method
	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	vtkPolyDataOrientNormalsByVoting(const vtkPolyDataOrientNormalsByVoting&);  // Not implemented.
	void operator=(const vtkPolyDataOrientNormalsByVoting&);  // Not implemented.

	//! Use a point locator to find the neighbouring points of all points
	void CreateConnectivity(vtkPolyData *input);

	//! Compute the cliques
	void ComputeCliques(vtkPolyData *input, vtkPolyData *output);
	
	//! Compute the cliques without using a pre-computed connectivity map
	void ComputeCliquesNoConnectivity(vtkPolyData *input, vtkPolyData *output);
		
	//! Check if two points are connected
	/** \param flip should normal be flipped */
	bool PointsConnected(int pid1, int pid2, vtkPolyData *input, vtkDataArray *normals, bool &flip);

	//! If a reference polydata is set, it is used to vote for the overall normal direction in each clique
	bool UseReferencePolyDataInVoting(vtkPolyData *input, vtkPolyData *output );

	//! Flip normals if the ends of the normals points somewhat inwards
	void MaxSpanningInVoting(vtkPolyData *input, vtkPolyData *output );
		
	//! Remove small cliques
	void CutCliques( vtkPolyData *input, vtkPolyData *output );

	//! Keep only largest clique
	void KeepLargestClique( vtkPolyData *input, vtkPolyData *output );

	//! Locator search radius for finding connectivity
	double SearchRadius;

	//! Keep track of neighbours of all points
	std::vector<std::vector<int> > m_Connectivity;

	//! The cliques of points
	std::vector<std::vector<int> > m_Cliques;

	//! Reference poly data
	vtkPolyData *ReferencePolyData;

	// The minimum clique that will be kept
	int MinCliqueSize;

	//! The max angle between two neighbouring points
	double NeighbourAngle;

	//! Only keep largest clique
	int LargestCliqueOnly;
};

#endif
