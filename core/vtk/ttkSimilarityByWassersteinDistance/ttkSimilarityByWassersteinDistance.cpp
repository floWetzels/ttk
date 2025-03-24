#include <ttkBottleneckDistanceUtils.h>
#include <ttkMacros.h>
#include <ttkSimilarityByWassersteinDistance.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkImageData.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByWassersteinDistance);

ttkSimilarityByWassersteinDistance::ttkSimilarityByWassersteinDistance() {
  this->setDebugMsgPrefix("ttkSimilarityByWassersteinDistance");
}

ttk::CriticalVertex* getExtremumVertex(ttk::PersistencePair& pair){
  if(pair.birth.type==ttk::CriticalType::Local_maximum || pair.birth.type==ttk::CriticalType::Local_minimum)
    return &pair.birth;
  else
    return &pair.death;
}

int ttkSimilarityByWassersteinDistance::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Number of input files
  auto input = vtkMultiBlockDataSet::GetData(inputVector[0]);

  const int numInputs = input->GetNumberOfBlocks();

  std::vector<vtkUnstructuredGrid *> inputVTUs(numInputs);
  for(int i = 0; i < numInputs; i++)
    inputVTUs[i] = vtkUnstructuredGrid::SafeDownCast(input->GetBlock(i));

  std::vector<ttk::DiagramType> inputPersistenceDiagrams(numInputs);
  std::vector<vtkNew<vtkUnstructuredGrid>> outputDiags(2 * numInputs - 2);
  std::vector<std::vector<ttk::MatchingType>> outputMatchings(numInputs - 1);

  std::string const algorithm = DistanceAlgorithm;
  double const tolerance = Tolerance;
  std::string const wasserstein = WassersteinMetric;

  for(int i = 0; i < numInputs; ++i)
    VTUToDiagram(inputPersistenceDiagrams[i], inputVTUs[i], *this);

  {
    ttk::Timer timer;
    this->printMsg("Perform Matchings", 0, 0, ttk::debug::LineMode::REPLACE);
    this->performMatchings(
      numInputs, inputPersistenceDiagrams, outputMatchings,
      algorithm, // Not from paraview, from enclosing tracking plugin
      wasserstein, tolerance, PX, PY, PZ, PS, PE // Coefficients
    );
    this->printMsg("Perform Matchings", 1, timer.getElapsedTime());
  }


  // Get back meshes.
  for(int i = 0; i < numInputs - 1; ++i) {
    outputDiags[2 * i + 0]->ShallowCopy(inputVTUs[i]);
    outputDiags[2 * i + 1]->ShallowCopy(inputVTUs[i + 1]);

    int const status = augmentDiagrams(
      outputMatchings[i], outputDiags[2 * i + 0], outputDiags[2 * i + 1]);
    if(status < 0)
      return -2;
  }

  for(int i = 0; i < numInputs; ++i) {
    const auto &grid = outputDiags[i];
    if(!grid || !grid->GetCellData()
      || !grid->GetCellData()->GetArray("Persistence")) {
      this->printErr("Inputs are not persistence diagrams");
      return 0;
    }

    // Check if inputs have the same data type and the same number of points
    if(grid->GetCellData()->GetArray("Persistence")->GetDataType()
      != outputDiags[0]
            ->GetCellData()
            ->GetArray("Persistence")
            ->GetDataType()) {
      this->printErr("Inputs of different data types");
      return 0;
    }
  }

  for(int i = 0; i < numInputs - 1; ++i) {
    const auto &grid1 = outputDiags[i];
    const auto &grid2 = outputDiags[i + 1];
    if(i % 2 == 1 && i < numInputs - 1
      && grid1->GetCellData()->GetNumberOfTuples()
            != grid2->GetCellData()->GetNumberOfTuples()) {
      this->printErr("Inconsistent length or order of input diagrams.");
      return 0;
    }
  }

  std::vector<ttk::trackingTuple> trackings;
  {
    ttk::Timer timer;
    this->printMsg("Perform Tracking", 0, 0, ttk::debug::LineMode::REPLACE);
    this->performTracking(inputPersistenceDiagrams, outputMatchings, trackings);
    this->printMsg("Perform Tracking", 1, timer.getElapsedTime());
  }

  auto prepArray = [](int size){
    auto array = vtkSmartPointer<vtkIntArray>::New();
    array->SetName("VertexId");
    array->SetNumberOfTuples(size);
    return array;
  };

  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  std::vector<std::unordered_map<int, int>> idIndexMaps(numInputs);
  {
    ttk::Timer timer;
    this->printMsg("Initialize Output", 0, 0, ttk::debug::LineMode::REPLACE);

    int nf0 = inputPersistenceDiagrams[0].size();
    for(int t=1; t<numInputs; t++) {
      int nf1 = inputPersistenceDiagrams[t].size();
      auto matrix = ttkSimilarityAlgorithm::InitializeMatrix(
        "Matching",

        VTK_INT,
        nf0, nf1
      );
      matrix->GetPointData()->GetArray(0)->Fill(0);

      auto indexIdMap0 = prepArray(nf0);
      auto indexIdMap1 = prepArray(nf1);
      for(int p=0; p<inputPersistenceDiagrams[t-1].size(); p++)
        indexIdMap0->SetValue(p,getExtremumVertex(inputPersistenceDiagrams[t-1][p])->id);
      for(int p=0; p<inputPersistenceDiagrams[t].size(); p++)
        indexIdMap1->SetValue(p,getExtremumVertex(inputPersistenceDiagrams[t][p])->id);

      if(t==1)
        ttkSimilarityAlgorithm::BuildIdIndexMap(idIndexMaps[0],indexIdMap0);
      ttkSimilarityAlgorithm::BuildIdIndexMap(idIndexMaps[t],indexIdMap1);

      ttkSimilarityAlgorithm::AddIndexIdMaps(matrix, indexIdMap0, indexIdMap1);

      output->SetBlock(t-1,matrix);
      nf0 = nf1;
    }
    this->printMsg("Initialize Output", 1, timer.getElapsedTime());
  }

  {
    ttk::Timer timer;
    this->printMsg("Initialize Output", 0, 0, ttk::debug::LineMode::REPLACE);

    for(size_t i=0; i<trackings.size(); i++){
      const ttk::trackingTuple &tt = trackings[i];
      const int numStart = std::get<0>(tt);
      const std::vector<ttk::SimplexId> &chain = std::get<2>(tt);
      int const chainLength = chain.size();
      // std::cout<<"-------------------------------"<<std::endl;
      // std::cout<<i<<" ["<<numStart<<","<<numEnd<<","<<chainLength<<"] ";
      for(int c=0; c<chainLength-1; c++){
        int const t0 = numStart + c;
        int const t1 = t0 + 1; // c % 2 == 0 ? d1id + 1 : d1id;

        auto matrixData = ttkUtils::GetPointer<int>(
          static_cast<vtkImageData*>(output->GetBlock(t0))->GetPointData()->GetArray(0)
        );
        auto& map0 = idIndexMaps[t0];
        auto& map1 = idIndexMaps[t1];
        auto &diagram1 = inputPersistenceDiagrams[t0];
        auto &diagram2 = inputPersistenceDiagrams[t1];

        const auto n1 = chain[c];
        const auto n2 = chain[c + 1];
        auto v0 = getExtremumVertex(diagram1[n1]);
        auto v1 = getExtremumVertex(diagram2[n2]);
        auto v0Idx = map0.find(v0->id);
        auto v1Idx = map1.find(v1->id);
        if(v0Idx!=map0.end() && v1Idx!=map1.end())
          matrixData[v1Idx->second*map0.size()+v0Idx->second] = 1;
          // matrixData[v0Idx->second*map1.size()+v1Idx->second] = 1;
      }

      // std::cout<<std::endl;
    }

    this->printMsg("Initialize Output", 1, timer.getElapsedTime());
  }




  // // Build mesh.
  // buildMesh(trackingsBase, outputMatchings, inputPersistenceDiagrams,
  //           useGeometricSpacing, spacing, DoPostProc, trackingTupleToMerged,
  //           points, persistenceDiagram, persistenceScalars, valueScalars,
  //           matchingIdScalars, lengthScalars, timeScalars, componentIds,
  //           pointTypeScalars, *this);

  // outputMesh->ShallowCopy(persistenceDiagram);

  return 1;
}
