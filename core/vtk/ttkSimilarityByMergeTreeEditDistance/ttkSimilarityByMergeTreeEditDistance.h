#pragma once

// VTK Module
#include <ttkSimilarityByMergeTreeEditDistanceModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityByMergeTreeEditDistance.h>

class TTKSIMILARITYBYMERGETREEEDITDISTANCE_EXPORT ttkSimilarityByMergeTreeEditDistance
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByMergeTreeEditDistance {

private:
  bool NormalizeMatrix{true};
  int backend = 0;
  int branchMetric = 0;

public:
  static ttkSimilarityByMergeTreeEditDistance *New();
  vtkTypeMacro(ttkSimilarityByMergeTreeEditDistance, ttkSimilarityAlgorithm);

  vtkSetMacro(NormalizeMatrix, bool);
  vtkGetMacro(NormalizeMatrix, bool);

  void SetBranchMetric(int b) {
    branchMetric = b;
    Modified();
  }
  vtkGetMacro(branchMetric, int);

  void SetBackend(int b) {
    backend = b;
    Modified();
  }
  vtkGetMacro(backend, int);

protected:
  ttkSimilarityByMergeTreeEditDistance();
  ~ttkSimilarityByMergeTreeEditDistance();

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;
};
