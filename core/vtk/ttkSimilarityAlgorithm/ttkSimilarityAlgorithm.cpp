#include <ttkSimilarityAlgorithm.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkStringArray.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityAlgorithm);

ttkSimilarityAlgorithm::ttkSimilarityAlgorithm() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkSimilarityAlgorithm::~ttkSimilarityAlgorithm() {
}

int ttkSimilarityAlgorithm::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port >= 0 && port < this->GetNumberOfInputPorts()) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkSimilarityAlgorithm::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

template <typename DT>
void BuildIdIndexMapDT(std::unordered_map<ttk::SimplexId, ttk::SimplexId> &map,
                       const int n,
                       const DT *indexIdMapData) {
  for(int i = 0; i < n; i++)
    map.emplace(
      std::make_pair(static_cast<ttk::SimplexId>(indexIdMapData[i]), i));
}

int ttkSimilarityAlgorithm::BuildIdIndexMap(
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> &map, const vtkDataArray *indexIdMap) {
  if(!indexIdMap)
    return 0;

  switch(indexIdMap->GetDataType()) {
    vtkTemplateMacro(
      BuildIdIndexMapDT<VTK_TT>(map, indexIdMap->GetNumberOfValues(),
                                ttkUtils::GetConstPointer<VTK_TT>(indexIdMap)));
  }

  return 1;
}

int ttkSimilarityAlgorithm::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto output = vtkImageData::GetData(outputVector);
  auto input = vtkDataObject::GetData(inputVector[0]);

  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> sequence;

  auto iterationInformation = vtkDoubleArray::SafeDownCast(input->GetFieldData()->GetArray("_ttk_IterationInfo"));
  if(iterationInformation->GetValue(0)>0){
    if(!this->ComputeSimilarityMatrix(output,this->PreviousInputs,input))
        return 0;
  } else {
    vtkNew<vtkStringArray> ttk_ignore;
    ttk_ignore->SetName("_ttk_IterationIgnore");
    output->GetFieldData()->AddArray(ttk_ignore);
  }

  this->PreviousInputs = vtkSmartPointer<vtkDataObject>::Take(input->NewInstance());
  this->PreviousInputs->DeepCopy(input);

  return 1;
}

std::string ttkSimilarityAlgorithm::GetIdArrayName(vtkFieldData *fieldData) {
  std::string result;

  int found = 0;
  for(int a = 0; a < fieldData->GetNumberOfArrays(); a++) {
    auto array = fieldData->GetAbstractArray(a);
    const auto name = std::string(array->GetName());
    const int n = name.size();

    if(name.substr(std::max(0,n - 4)).compare("_t-1") == 0) {
      result = name.substr(0, n - 4);
      found++;
    } else if(name.substr(std::max(0,n - 2)).compare("_t") == 0) {
      found++;
    }
  }

  return found == 2 ? result : "";
};

int ttkSimilarityAlgorithm::GetIndexIdMaps(vtkDataArray *&indexIdMapP,
                                           vtkDataArray *&indexIdMapC,
                                           vtkFieldData *fieldData) {
  std::string idArrayName = ttkSimilarityAlgorithm::GetIdArrayName(fieldData);

  if(idArrayName.size() < 1)
    return 0;

  indexIdMapP = fieldData->GetArray((idArrayName + "_t-1").data());
  indexIdMapC = fieldData->GetArray((idArrayName + "_t").data());

  return 1;
};

int ttkSimilarityAlgorithm::AddIndexIdMap(vtkImageData *similarityMatrix,
                                          vtkDataArray *idArray,
                                          const bool isMapForCurrentTimestep) {
  auto fd = similarityMatrix->GetFieldData();
  auto array = vtkSmartPointer<vtkDataArray>::Take(idArray->NewInstance());
  array->ShallowCopy(idArray);
  array->SetName((std::string(idArray->GetName())
                  + (isMapForCurrentTimestep ? "_t" : "_t-1"))
                   .data());
  fd->AddArray(array);
  return 1;
}

int ttkSimilarityAlgorithm::AddIndexIdMaps(vtkImageData *similarityMatrix,
                                           vtkDataArray *indexIdMapP,
                                           vtkDataArray *indexIdMapC) {
  int status = 0;
  status = ttkSimilarityAlgorithm::AddIndexIdMap(
    similarityMatrix, indexIdMapP, false);
  if(!status)
    return 0;

  status = ttkSimilarityAlgorithm::AddIndexIdMap(
    similarityMatrix, indexIdMapC, true);
  if(!status)
    return 0;

  return 1;
}

int ttkSimilarityAlgorithm::AddIndexIdMaps(
  vtkImageData *similarityMatrix,
  const std::unordered_map<ttk::SimplexId, ttk::SimplexId> &idIndexMapP,
  const std::unordered_map<ttk::SimplexId, ttk::SimplexId> &idIndexMapC,
  const std::string &idArrayName) {
  auto fd = similarityMatrix->GetFieldData();
  int a = 0;
  const std::string suffix[2]{"_t-1", "_t"};
  for(auto map :
      std::vector<const std::unordered_map<ttk::SimplexId, ttk::SimplexId> *>(
        {&idIndexMapP, &idIndexMapC})) {
    auto array = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    array->SetName((idArrayName + suffix[a++]).data());
    array->SetNumberOfTuples(map->size());
    auto arrayData = ttkUtils::GetPointer<ttk::SimplexId>(array);
    for(auto it : (*map))
      arrayData[it.second] = it.first;

    fd->AddArray(array);
  }

  return 1;
}
