#include <ttkComponentSize.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkConnectivityFilter.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

vtkStandardNewMacro(ttkComponentSize);

ttkComponentSize::ttkComponentSize() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkComponentSize::~ttkComponentSize() = default;

int ttkComponentSize::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkComponentSize::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkComponentSize::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  const size_t nPoints = input->GetNumberOfPoints();

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  auto sIdArray = this->GetInputArrayToProcess(0, inputVector);
  auto sIdArrayData = ttkUtils::GetPointer<const int>(sIdArray);

  sIdArray->Print(std::cout);

  std::unordered_map<int,ttk::ComponentSize::Component> components;
  if(!this->initializeComponents(
    components,
    nPoints,
    sIdArrayData
  ))
    return 0;

  if(!this->computeComponents(
    components,
    nPoints,
    sIdArrayData,
    triangulation
  ))
    return 0;

  {
    ttk::Timer t;
    this->printMsg("Generating Output", 0, 0, ttk::debug::LineMode::REPLACE);

    const size_t nComponents = components.size();
    auto output = vtkPolyData::GetData(outputVector);
    // points
    {
      auto cSizeArray = vtkSmartPointer<vtkFloatArray>::New();
      cSizeArray->SetName("Size");
      cSizeArray->SetNumberOfTuples(nComponents);
      auto cSizeArrayData = ttkUtils::GetPointer<float>(cSizeArray);

      auto cIdArray = vtkSmartPointer<vtkIntArray>::New();
      cIdArray->SetName(sIdArray->GetName());
      cIdArray->SetNumberOfTuples(nComponents);
      auto cIdArrayData = ttkUtils::GetPointer<int>(cIdArray);

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetDataTypeToFloat();
      points->SetNumberOfPoints(nComponents);
      auto pointsData = ttkUtils::GetPointer<float>(points->GetData());

      int i=0;
      int j=0;
      for(const auto& it: components){
        pointsData[j++] = it.second.center[0];
        pointsData[j++] = it.second.center[1];
        pointsData[j++] = it.second.center[2];

        cIdArrayData[i] = it.first;
        cSizeArrayData[i] = it.second.size;

        i++;
      }
      output->SetPoints(points);

      auto pd = output->GetPointData();
      pd->AddArray(cIdArray);
      pd->AddArray(cSizeArray);
    }

    // cells
    {
      auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
      connectivityArray->SetNumberOfTuples(nComponents);
      auto connectivityArrayData = ttkUtils::GetPointer<int>(connectivityArray);
      for(size_t i = 0; i < nComponents; i++)
        connectivityArrayData[i] = i;

      auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
      offsetArray->SetNumberOfTuples(nComponents + 1);
      auto offsetArrayData = ttkUtils::GetPointer<int>(offsetArray);
      for(size_t i = 0; i <= nComponents; i++)
        offsetArrayData[i] = i;

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);

      output->SetVerts(cellArray);
    }

    // Copy Field Data
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    this->printMsg("Generating Output", 1, t.getElapsedTime(), this->getThreadNumber());
  }

  return 1;
}
