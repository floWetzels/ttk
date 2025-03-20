#include <ttkSimilarityByMergeTreeSegmentation.h>

#include <vtkObjectFactory.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByMergeTreeSegmentation);

ttkSimilarityByMergeTreeSegmentation::ttkSimilarityByMergeTreeSegmentation() {
  this->SetNumberOfInputPorts(2);
}

ttkSimilarityByMergeTreeSegmentation::~ttkSimilarityByMergeTreeSegmentation() {}

int ttkSimilarityByMergeTreeSegmentation::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto inputMergeTrees = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto inputSegmentations = vtkMultiBlockDataSet::GetData(inputVector[1]);
  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  const size_t n = inputMergeTrees->GetNumberOfBlocks();

  for(size_t t=1; t<n; t++){
    auto m0 = vtkDataSet::SafeDownCast(inputMergeTrees->GetBlock(t-1));
    auto m1 = vtkDataSet::SafeDownCast(inputMergeTrees->GetBlock(t));

    auto d0 = vtkDataSet::SafeDownCast(inputSegmentations->GetBlock(t-1));
    auto d1 = vtkDataSet::SafeDownCast(inputSegmentations->GetBlock(t));

    const int nNodes0 = m0->GetNumberOfPoints();
    const int nNodes1 = m1->GetNumberOfPoints();

    // extract arrays
    auto seg0 = vtkIntArray::SafeDownCast(d0->GetPointData()->GetArray("NodeId"));
    auto seg1 = vtkIntArray::SafeDownCast(d1->GetPointData()->GetArray("NodeId"));
    if(!seg0 || !seg1)
      return !this->printErr(
        "Unable to retrieve `NodeId` arrays from segmentations.");

    auto next0
      = vtkIntArray::SafeDownCast(m0->GetPointData()->GetArray("NextId"));
    auto next1
      = vtkIntArray::SafeDownCast(m1->GetPointData()->GetArray("NextId"));
    if(!next0 || !next1)
      return !this->printErr(
        "Unable to retrieve `NextId` arrays from merge trees.");

    auto scalars0 = this->GetInputArrayToProcess(0, m0);
    auto scalars1 = this->GetInputArrayToProcess(0, m1);
    if(!scalars0 || !scalars1)
      return !this->printErr("Unable to retrieve merge tree scalar arrays.");

    // initialize similarity matrix
    auto matrix = ttkSimilarityAlgorithm::InitializeMatrix(
      "Overlap", VTK_INT, nNodes0, nNodes1
    );
    auto matrixData = matrix->GetPointData()->GetArray(0);

    // compute overlap of segments
    int status = 0;
    ttkTypeMacroA(
      scalars0->GetDataType(),
      (status = this->computeSegmentationOverlap<int, T0>(
        ttkUtils::GetPointer<int>(matrixData),

        ttkUtils::GetConstPointer<int>(seg0),
        ttkUtils::GetConstPointer<int>(seg1), seg0->GetNumberOfTuples(),
        ttkUtils::GetConstPointer<int>(next0),
        ttkUtils::GetConstPointer<int>(next1),
        ttkUtils::GetConstPointer<T0>(scalars0),
        ttkUtils::GetConstPointer<T0>(scalars1), nNodes0, nNodes1)));
    if(!status)
      return 0;

    status = ttkSimilarityAlgorithm::AddIndexIdMaps(
      matrix, this->GetInputArrayToProcess(1, m0),
      this->GetInputArrayToProcess(1, m1));
    if(!status)
      return 0;

    output->SetBlock(t-1,matrix);
  }

  return 1;
}
