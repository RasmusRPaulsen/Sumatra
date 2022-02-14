#ifndef __vtkPolyDataProjection_h
#define __vtkPolyDataProjection_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkCellLocator.h"

// Projects the source surface on the target surface

class vtkPolyDataProjection : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPolyDataProjection,vtkPolyDataAlgorithm);
  static vtkPolyDataProjection *New();
  void PrintSelf(ostream& os, vtkIndent indent);

	void SetTargetData(vtkPolyData *pd);
	vtkPolyData *GetTarget();

	vtkPolyData *GetDisplacementField();

	vtkSetMacro(CreateDisplacementField, bool);
	vtkGetMacro(CreateDisplacementField, bool);

protected:
	vtkPolyDataProjection();
	~vtkPolyDataProjection();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
   virtual int FillInputPortInformation(int port, vtkInformation *info);

	bool CreateDisplacementField;

	vtkPolyData *m_DisplacementField;

private:
	vtkPolyDataProjection(const vtkPolyDataProjection&) {};
	void operator=(const vtkPolyDataProjection&) {};
};

#endif
