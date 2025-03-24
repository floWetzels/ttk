#pragma once

// VTK Module
#include <ttkSimilarityByJacobiSetModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

class TTKSIMILARITYBYJACOBISET_EXPORT ttkSimilarityByJacobiSet
  : public ttkSimilarityAlgorithm {

public:
  static ttkSimilarityByJacobiSet *New();
  vtkTypeMacro(ttkSimilarityByJacobiSet, ttkSimilarityAlgorithm);

protected:
  ttkSimilarityByJacobiSet();
  ~ttkSimilarityByJacobiSet();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
