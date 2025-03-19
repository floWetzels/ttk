#include <ttkSimilarityByMergeTreeEditDistance.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByMergeTreeEditDistance);

ttkSimilarityByMergeTreeEditDistance::ttkSimilarityByMergeTreeEditDistance() {}
ttkSimilarityByMergeTreeEditDistance::~ttkSimilarityByMergeTreeEditDistance() {}

int ttkSimilarityByMergeTreeEditDistance::ComputeSimilarityMatrix(
  vtkImageData *similarityMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {

  // unpack input
  auto t0 = vtkMultiBlockDataSet::SafeDownCast(inputDataObjects0);
  auto t1 = vtkMultiBlockDataSet::SafeDownCast(inputDataObjects1);
  if(!t0 || !t1)
    return !this->printErr("Input data objects need to be vtkPointSets.");
  
  vtkNew<vtkMultiBlockDataSet> mtmb;
  vtkNew<vtkMultiBlockDataSet> mtnodes;
  vtkNew<vtkMultiBlockDataSet> mtarcs;

  mtnodes->SetNumberOfBlocks(2);
  mtnodes->SetBlock(0, t0->GetBlock(0));
  mtnodes->SetBlock(1, t1->GetBlock(0));

  mtarcs->SetNumberOfBlocks(2);
  mtarcs->SetBlock(0, t0->GetBlock(1));
  mtarcs->SetBlock(1, t1->GetBlock(1));

  mtmb->SetNumberOfBlocks(2);
  mtmb->SetBlock(0, mtnodes);
  mtmb->SetBlock(1, mtarcs);

  auto nodes0 = vtkUnstructuredGrid::SafeDownCast(t0->GetBlock(0));
  auto nodes1 = vtkUnstructuredGrid::SafeDownCast(t1->GetBlock(0));

  vtkNew<ttkMergeTreeFeatureTracking> ft;
  ft->SetBackend(backend);
  if(this->backend==3){
    ft->SetBranchMetric(branchMetric);
  }

  ft->SetInputDataObject(0,mtmb.GetPointer());
  ft->SetImportantPairs(0);
  if(this->backend==0){
    ft->SetEpsilonTree1(5);
  } else {
    ft->SetEpsilonTree1(0);
  }
  ft->SetEpsilon2Tree1(100);
  ft->SetEpsilon3Tree1(100);
  //ft->SetPlanarLayout(true);
  //ft->SetNormalizedWasserstein(false);
  ft->Update();

  // get point coordinates
  auto coords0 = nodes0->GetPoints()->GetData();
  auto coords1 = nodes1->GetPoints()->GetData();

  // initialize similarity matrix i.e., distance matrix
  similarityMatrix->SetDimensions(nodes0->GetNumberOfPoints(), nodes1->GetNumberOfPoints(), 1);
  // similarityMatrix->AllocateScalars(coords0->GetDataType(), 1);
  auto matrixData = similarityMatrix->GetPointData()->GetArray(0);
  matrixData->SetName("Matched");
  matrixData->SetNumberOfTuples(nodes0->GetNumberOfPoints()*nodes1->GetNumberOfPoints());
  matrixData->SetNumberOfComponents(1);

  auto matchings_mb = vtkMultiBlockDataSet::New();
  matchings_mb->ShallowCopy(vtkMultiBlockDataSet::SafeDownCast(ft->GetOutputDataObject(1)));
  for(ttk::SimplexId i=0; i<matchings_mb->GetNumberOfBlocks(); i++){
    auto matchingi_vtk = vtkUnstructuredGrid::SafeDownCast(matchings_mb->GetBlock(i));
    matchings[i] = std::vector<ttk::SimplexId>(currmemberNodes->GetNumberOfPoints(),-1);
    for(ttk::SimplexId cellIdx = 0; cellIdx < matchingi_vtk->GetNumberOfCells(); cellIdx++) {
      ttk::SimplexId id1 = matchingi_vtk->GetCellData()->GetArray("tree1NodeId")->GetComponent(cellIdx,0);
      ttk::SimplexId id2 = matchingi_vtk->GetCellData()->GetArray("tree2NodeId")->GetComponent(cellIdx,0);
      ttk::SimplexId n1 = nodes1->GetPointData()->GetArray("NodeId")->GetComponent(id1,0);
      ttk::SimplexId n2 = nodes2->GetPointData()->GetArray("NodeId")->GetComponent(id2,0);
      matrixData->SetTuple1(id1 + id2*nodes0->GetNumberOfPoints());
    }
  }

  return 1;
}
