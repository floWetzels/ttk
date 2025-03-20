#pragma once

// VTK Module
#include <ttkSimilarityByDistanceModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityByDistance.h>

class TTKSIMILARITYBYDISTANCE_EXPORT ttkSimilarityByDistance
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByDistance {

private:
  bool NormalizeMatrix{true};

public:
  static ttkSimilarityByDistance *New();
  vtkTypeMacro(ttkSimilarityByDistance, ttkSimilarityAlgorithm);

  vtkSetMacro(NormalizeMatrix, bool);
  vtkGetMacro(NormalizeMatrix, bool);

protected:
  ttkSimilarityByDistance();
  ~ttkSimilarityByDistance();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

};
