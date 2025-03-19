#pragma once

// VTK Module
#include <ttkMergeTreeRemakeModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <MergeTreeRemake.h>

class TTKMERGETREEREMAKE_EXPORT ttkMergeTreeRemake : public ttkAlgorithm,
                                         public ttk::mt::MergeTreeRemake {
private:
  int Type{0};

public:
  vtkGetMacro(Type, int);
  vtkSetMacro(Type, int);

  static ttkMergeTreeRemake *New();
  vtkTypeMacro(ttkMergeTreeRemake, ttkAlgorithm);

protected:
  ttkMergeTreeRemake();
  ~ttkMergeTreeRemake() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};