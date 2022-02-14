#ifndef __vtkPolyDataTextReader_h
#define __vtkPolyDataTextReader_h

#include "vtkPolyDataAlgorithm.h"

//! Reader for text files
/** Can read points and normals. Tries to automatically determine delimiter (default is space) */
class vtkPolyDataTextReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkPolyDataTextReader *New();
  vtkTypeMacro(vtkPolyDataTextReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of  file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  bool ReaderError() const;

protected:
  vtkPolyDataTextReader();
  ~vtkPolyDataTextReader();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
private:
  vtkPolyDataTextReader(const vtkPolyDataTextReader&);  // Not implemented.
  void operator=(const vtkPolyDataTextReader&);  // Not implemented.

  bool ErrorReading;
};

#endif
