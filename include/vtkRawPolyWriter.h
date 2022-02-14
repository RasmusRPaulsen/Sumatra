#ifndef __vtkRawPolyWriter_h
#define __vtkRawPolyWriter_h

#include "vtkPolyDataWriter.h"

class vtkRawPolyWriter : public vtkPolyDataWriter
{
public:
  static vtkRawPolyWriter *New();
  vtkTypeMacro(vtkRawPolyWriter,vtkPolyDataWriter);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Write scalars if they are present
  vtkSetMacro(WriteScalars, int);
  vtkGetMacro(WriteScalars, int);
  vtkBooleanMacro(WriteScalars, int);

  // Description:
  // Write normals if they are present
  vtkSetMacro(WriteNormals, int);
  vtkGetMacro(WriteNormals, int);
  vtkBooleanMacro(WriteNormals, int);

protected:
  vtkRawPolyWriter();
  ~vtkRawPolyWriter() {};

  void WriteData();

  void WriteAsciiRAW(vtkPolyData *pd);


private:
  vtkRawPolyWriter(const vtkRawPolyWriter&);  // Not implemented.
  void operator=(const vtkRawPolyWriter&);  // Not implemented.

  int WriteScalars;

  int WriteNormals;

};

#endif

