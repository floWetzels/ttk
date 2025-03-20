#include <ttkSimilarityByMergeTreeEditDistance.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <ttkMergeTreeFeatureTracking.h>

vtkStandardNewMacro(ttkSimilarityByMergeTreeEditDistance);

ttkSimilarityByMergeTreeEditDistance::ttkSimilarityByMergeTreeEditDistance() {
  this->SetNumberOfInputPorts(4);
}
ttkSimilarityByMergeTreeEditDistance::~ttkSimilarityByMergeTreeEditDistance() {}

int ttkSimilarityByMergeTreeEditDistance::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto inputNodes = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto inputArcs = vtkMultiBlockDataSet::GetData(inputVector[1]);
  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  const size_t n = inputNodes->GetNumberOfBlocks();

  for(size_t t=1; t<n; t++){
    vtkNew<vtkMultiBlockDataSet> mtmb;
    vtkNew<vtkMultiBlockDataSet> mtnodes;
    vtkNew<vtkMultiBlockDataSet> mtarcs;

    mtnodes->SetNumberOfBlocks(2);
    mtnodes->SetBlock(0, inputNodes->GetBlock(t-1));
    mtnodes->SetBlock(1, inputNodes->GetBlock(t));

    mtarcs->SetNumberOfBlocks(2);
    mtarcs->SetBlock(0, inputArcs->GetBlock(t-1));
    mtarcs->SetBlock(1, inputArcs->GetBlock(t));

    mtmb->SetNumberOfBlocks(2);
    mtmb->SetBlock(0, mtnodes);
    mtmb->SetBlock(1, mtarcs);

    auto nodes0 = vtkUnstructuredGrid::SafeDownCast(inputNodes->GetBlock(t-1));
    auto nodes1 = vtkUnstructuredGrid::SafeDownCast(inputNodes->GetBlock(t));

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

    // initialize similarity matrix i.e., distance matrix
    auto matrix = ttkSimilarityAlgorithm::InitializeMatrix(
      "Matched", VTK_INT, nodes0->GetNumberOfPoints(), nodes1->GetNumberOfPoints()
    );
    auto matrixData = matrix->GetPointData()->GetArray(0);
    for(ttk::SimplexId i=0; i<nodes0->GetNumberOfPoints()*nodes1->GetNumberOfPoints(); i++)
      matrixData->SetTuple1(i,0);

    vtkNew<vtkMultiBlockDataSet> matchings_mb;
    matchings_mb->ShallowCopy(vtkMultiBlockDataSet::SafeDownCast(ft->GetOutputDataObject(1)));
    for(ttk::SimplexId i=0; i<matchings_mb->GetNumberOfBlocks(); i++){
      auto matchingi_vtk = vtkUnstructuredGrid::SafeDownCast(matchings_mb->GetBlock(i));
      for(ttk::SimplexId cellIdx = 0; cellIdx < matchingi_vtk->GetNumberOfCells(); cellIdx++) {
        ttk::SimplexId id1 = matchingi_vtk->GetCellData()->GetArray("tree1NodeId")->GetComponent(cellIdx,0);
        ttk::SimplexId id2 = matchingi_vtk->GetCellData()->GetArray("tree2NodeId")->GetComponent(cellIdx,0);
        ttk::SimplexId n1 = nodes0->GetPointData()->GetArray("NodeId")->GetComponent(id1,0);
        ttk::SimplexId n2 = nodes1->GetPointData()->GetArray("NodeId")->GetComponent(id2,0);
        matrixData->SetTuple1(id1 + id2*nodes0->GetNumberOfPoints(),1);
      }
    }
    int status = 0;

    auto indexIdMap0 = this->GetInputArrayToProcess(0, nodes0);
    auto indexIdMap1 = this->GetInputArrayToProcess(0, nodes1);
    if(!indexIdMap0 || !indexIdMap1)
      return !this->printErr("Unable to retrieve feature IDs.");

    status = ttkSimilarityAlgorithm::AddIndexIdMaps(
      matrix, indexIdMap0, indexIdMap1);
    if(!status)
      return 0;

    output->SetBlock(t-1,matrix);
  }

  return 1;
}
