#include <ttkSimilarityByOverlap.h>

#include <vtkObjectFactory.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByOverlap);
ttkSimilarityByOverlap::ttkSimilarityByOverlap() {}
ttkSimilarityByOverlap::~ttkSimilarityByOverlap() {}

int ttkSimilarityByOverlap::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto input = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  const size_t n = input->GetNumberOfBlocks();

  for(size_t t=1; t<n; t++){
    // get id arrays
    auto ids0 = this->GetInputArrayToProcess(0, input->GetBlock(t-1));
    auto ids1 = this->GetInputArrayToProcess(0, input->GetBlock(t));

    // validate arrays
    if(!ids0 || !ids1)
      return !this->printErr("Unable to retrieve ids.");

    if(ids0->GetNumberOfComponents() != 1 || ids1->GetNumberOfComponents() != 1)
      return !this->printErr("Ids must have exactly one component.");

    if(ids0->GetNumberOfTuples() != ids1->GetNumberOfTuples())
      return !this->printErr("Ids must have same number of values.");

    if(ids0->GetDataType() != ids1->GetDataType())
      return !this->printErr("Ids must have same data type.");

    const int nVertices = ids0->GetNumberOfTuples();

    // extract unique ids from volume
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> idIndexMap0;
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> idIndexMap1;
    int status = 0;
    for(auto &it : std::vector<std::pair<
          vtkDataArray *, std::unordered_map<ttk::SimplexId, ttk::SimplexId> *>>(
          {{ids0, &idIndexMap0}, {ids1, &idIndexMap1}})) {
      ttkTypeMacroA(
        ids0->GetDataType(),
        (status = this->computeIdIndexMap<T0, ttk::SimplexId>(
           *it.second, ttkUtils::GetPointer<const T0>(it.first), nVertices)));
      if(!status)
        return 0;
    }

    const int nIds0 = idIndexMap0.size();
    const int nIds1 = idIndexMap1.size();

    // initialize similarity matrix
    auto matrix = ttkSimilarityAlgorithm::InitializeMatrix(
      "Overlap", VTK_INT, nIds0, nIds1
    );
    auto matrixData = matrix->GetPointData()->GetArray(0);

    // compute overlaps
    ttkTypeMacroA(ids0->GetDataType(),
                  (status = this->computeAdjacencyMatrix<T0, ttk::SimplexId>(
                     ttkUtils::GetPointer<int>(matrixData),
                     ttkUtils::GetPointer<const T0>(ids0),
                     ttkUtils::GetPointer<const T0>(ids1), nVertices, idIndexMap0,
                     idIndexMap1)));
    if(!status)
      return 0;

    status = ttkSimilarityAlgorithm::AddIndexIdMaps(
      matrix, idIndexMap0, idIndexMap1, ids0->GetName());
    if(!status)
      return 0;

    output->SetBlock(t-1,matrix);
  }

  return 1;
};
