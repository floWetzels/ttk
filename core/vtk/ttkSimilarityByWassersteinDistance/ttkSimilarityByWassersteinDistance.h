#pragma once

#include <TrackingFromPersistenceDiagrams.h>
#include <ttkSimilarityAlgorithm.h>

// VTK Module
#include <ttkSimilarityByWassersteinDistanceModule.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

class vtkUnstructuredGrid;

class TTKSIMILARITYBYWASSERSTEINDISTANCE_EXPORT
  ttkSimilarityByWassersteinDistance
  : public ttkSimilarityAlgorithm,
    protected ttk::TrackingFromPersistenceDiagrams {

public:
  static ttkSimilarityByWassersteinDistance *New();

  vtkTypeMacro(ttkSimilarityByWassersteinDistance, ttkSimilarityAlgorithm);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, const std::string &);
  vtkGetMacro(DistanceAlgorithm, std::string);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  // static int buildMesh(
  //   const std::vector<ttk::trackingTuple> &trackings,
  //   const std::vector<std::vector<ttk::MatchingType>> &outputMatchings,
  //   const std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
  //   const bool useGeometricSpacing,
  //   const double spacing,
  //   const bool doPostProc,
  //   const std::vector<std::set<int>> &trackingTupleToMerged,
  //   vtkPoints *points,
  //   vtkUnstructuredGrid *persistenceDiagram,
  //   vtkDoubleArray *persistenceScalars,
  //   vtkDoubleArray *valueScalars,
  //   vtkIntArray *matchingIdScalars,
  //   vtkIntArray *lengthScalars,
  //   vtkIntArray *timeScalars,
  //   vtkIntArray *componentIds,
  //   vtkIntArray *pointTypeScalars,
  //   const ttk::Debug &dbg);

protected:
  ttkSimilarityByWassersteinDistance();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // Input bottleneck config.
  double PostProcThresh{0.0};
  double Tolerance{1.0};
  double PX{1};
  double PY{1};
  double PZ{1};
  double PE{1};
  double PS{1};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"1"};
};
