/// \ingroup base
/// \class ttk::MergeTreeDistanceMatrix
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \author Florian Wetzels (wetzels@cs.uni-kl.de)
/// \date 2021.
///
/// This VTK filter uses the ttk::MergeTreeDistanceMatrix module to compute the
/// distance matrix of a group of merge trees.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \b Related \b publication \n
/// "Edit Distance between Merge Trees" \n
/// R. Sridharamurthy, T. B. Masood, A. Kamakshidasan and V. Natarajan. \n
/// IEEE Transactions on Visualization and Computer Graphics, 2020.
///
/// \b Related \b publication \n
/// "Branch Decomposition-Independent Edit Distances for Merge Trees." \n
/// Florian Wetzels, Heike Leitte, and Christoph Garth. \n
/// Computer Graphics Forum, 2022.
///
/// \b Related \b publication \n
/// "A Deformation-based Edit Distance for Merge Trees" \n
/// Florian Wetzels, Christoph Garth. \n
/// TopoInVis 2022.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// ttk common includes
#include <Debug.h>

#include <BranchMappingDistance.h>
#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <MergeTreeBase.h>
#include <MergeTreeDistance.h>
#include <PathMappingDistance.h>
#include <NaiveOneDegMergeTreeEditDistance.h>
#include <NaiveConstrMergeTreeEditDistance.h>

namespace ttk {

  /**
   * The MergeTreeDistanceMatrix class provides methods to compute the
   * distance between multiple merge trees and output a distance matrix.
   */
  class MergeTreeDistanceMatrix : virtual public Debug,
                                  virtual public MergeTreeBase {
  protected:
    int baseModule_ = 0;
    int branchMetric_ = 0;
    int pathMetric_ = 0;
    int pathMappingLookahead_ = 0;
    std::vector<std::map<ftm::idNode,int>> vid_to_seg;
    std::vector<std::map<double,int>> val_to_seg;

  public:
    MergeTreeDistanceMatrix() {
      this->setDebugMsgPrefix(
        "MergeTreeDistanceMatrix"); // inherited from Debug: prefix will be
                                    // printed at the
      // beginning of every msg
    }
    ~MergeTreeDistanceMatrix() override = default;

    void setBaseModule(int m) {
      baseModule_ = m;
    }

    void setBranchMetric(int m) {
      branchMetric_ = m;
    }

    void setPathMetric(int m) {
      pathMetric_ = m;
    }

    void setPathMappingLookahead(int l) {
      pathMappingLookahead_ = l;
    }

    /**
     * Implementation of the algorithm.
     */
    template <class dataType>
    void execute(std::vector<ftm::MergeTree<dataType>> &trees,
                 std::vector<ftm::MergeTree<dataType>> &trees2,
                 std::vector<std::vector<double>> &distanceMatrix) {
      treesNodeCorr_.resize(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        preprocessingPipeline<dataType>(
          trees[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          baseModule_ == 0 ? branchDecomposition_ : false, useMinMaxPair_, true,
          treesNodeCorr_[i], true, baseModule_ == 2);
      }
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      executePara<dataType>(trees, distanceMatrix);
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
      std::cout << "Matrix time: " << (float)time/1000. << "s" << std::endl;
      if(trees2.size() != 0) {
        std::vector<std::vector<int>> trees2NodeCorr(trees2.size());
        for(unsigned int i = 0; i < trees.size(); ++i) {
          preprocessingPipeline<dataType>(
            trees2[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
            baseModule_ == 0 ? branchDecomposition_ : false, useMinMaxPair_,
            true, treesNodeCorr_[i], true, baseModule_ == 2);
        }
        useDoubleInput_ = true;
        std::vector<std::vector<double>> distanceMatrix2(
          trees2.size(), std::vector<double>(trees2.size()));
        executePara<dataType>(trees2, distanceMatrix2, false);
        mixDistancesMatrix(distanceMatrix, distanceMatrix2);
      }
    }

    template <class dataType>
    void executePara(std::vector<ftm::MergeTree<dataType>> &trees,
                     std::vector<std::vector<double>> &distanceMatrix,
                     bool isFirstInput = true) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp single nowait
#endif
        executeParaImpl<dataType>(trees, distanceMatrix, isFirstInput);
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void executeParaImpl(std::vector<ftm::MergeTree<dataType>> &trees,
                         std::vector<std::vector<double>> &distanceMatrix,
                         bool isFirstInput = true) {
      std::cout << "------------\n";
      std::srand(17);
      int avg_time = 0;
      int avg_size = 0;
      unsigned int nruns = 100;
      std::vector<int> r1s;
      std::vector<int> r2s;
      for(unsigned int i = 0; i < nruns; ++i){
        auto r1 = std::rand() % distanceMatrix.size();
        auto r2 = std::rand() % distanceMatrix.size();
        r1s.push_back(r1);
        r2s.push_back(r2);
        int ttime = -1;
        float dist = -1;
        // std::cout << r1 << " " << r2 << " ; "
        //           << trees[r1].tree.getNumberOfNodes() << " "
        //           << trees[r2].tree.getNumberOfNodes()  << " ; ";

        if(baseModule_ == 0) {
          MergeTreeDistance mergeTreeDistance;
          mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
          mergeTreeDistance.setEpsilonTree1(epsilonTree1_);
          mergeTreeDistance.setEpsilonTree2(epsilonTree2_);
          mergeTreeDistance.setEpsilon2Tree1(epsilon2Tree1_);
          mergeTreeDistance.setEpsilon2Tree2(epsilon2Tree2_);
          mergeTreeDistance.setEpsilon3Tree1(epsilon3Tree1_);
          mergeTreeDistance.setEpsilon3Tree2(epsilon3Tree2_);
          mergeTreeDistance.setBranchDecomposition(branchDecomposition_);
          mergeTreeDistance.setParallelize(false);
          mergeTreeDistance.setPersistenceThreshold(persistenceThreshold_);
          mergeTreeDistance.setDebugLevel(std::min(debugLevel_, 2));
          mergeTreeDistance.setThreadNumber(this->threadNumber_);
          mergeTreeDistance.setNormalizedWasserstein(
            normalizedWasserstein_);
          mergeTreeDistance.setKeepSubtree(keepSubtree_);
          mergeTreeDistance.setDistanceSquaredRoot(distanceSquaredRoot_);
          mergeTreeDistance.setUseMinMaxPair(useMinMaxPair_);
          mergeTreeDistance.setPreprocess(false);
          // mergeTreeDistance.setSaveTree(true);
          mergeTreeDistance.setSaveTree(false);
          mergeTreeDistance.setCleanTree(true);
          mergeTreeDistance.setIsCalled(true);
          mergeTreeDistance.setPostprocess(false);
          mergeTreeDistance.setIsPersistenceDiagram(isPersistenceDiagram_);
          if(useDoubleInput_) {
            double const weight
              = mixDistancesMinMaxPairWeight(isFirstInput);
            mergeTreeDistance.setMinMaxPairWeight(weight);
            mergeTreeDistance.setDistanceSquaredRoot(true);
          }
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> outputMatching;
          std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
          dist = mergeTreeDistance.execute<dataType>(
            trees[r1], trees[r2], outputMatching);
          std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
          ttime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        } else if(baseModule_ == 1) {
          BranchMappingDistance branchDist;
          branchDist.setBaseMetric(branchMetric_);
          branchDist.setAssignmentSolver(assignmentSolverID_);
          branchDist.setSquared(distanceSquaredRoot_);
          branchDist.setEpsilonTree1(epsilonTree1_);
          branchDist.setEpsilonTree2(epsilonTree2_);
          branchDist.setEpsilon2Tree1(epsilon2Tree1_);
          branchDist.setEpsilon2Tree2(epsilon2Tree2_);
          branchDist.setEpsilon3Tree1(epsilon3Tree1_);
          branchDist.setEpsilon3Tree2(epsilon3Tree2_);
          branchDist.setPersistenceThreshold(persistenceThreshold_);
          branchDist.setPreprocess(false);
          // branchDist.setSaveTree(true);
          branchDist.setSaveTree(false);
          std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
          dist = branchDist.execute<dataType>(trees[r1], trees[r2]);
          std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
          ttime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        } else if(baseModule_ == 2) {
          PathMappingDistance pathDist;
          pathDist.setBaseMetric(pathMetric_);
          pathDist.setAssignmentSolver(assignmentSolverID_);
          pathDist.setSquared(distanceSquaredRoot_);
          pathDist.setComputeMapping(true);
          pathDist.setEpsilonTree1(epsilonTree1_);
          pathDist.setEpsilonTree2(epsilonTree2_);
          pathDist.setEpsilon2Tree1(epsilon2Tree1_);
          pathDist.setEpsilon2Tree2(epsilon2Tree2_);
          pathDist.setEpsilon3Tree1(epsilon3Tree1_);
          pathDist.setEpsilon3Tree2(epsilon3Tree2_);
          pathDist.setPersistenceThreshold(persistenceThreshold_);
          pathDist.setPreprocess(false);
          // pathDist.setSaveTree(true);
          pathDist.setSaveTree(false);
          pathDist.setlookahead(pathMappingLookahead_);
          std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
          dist = pathDist.execute<dataType>(trees[r1], trees[r2]);
          std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
          ttime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
          // std::cout << i << " " << j << " ; "
          //           << trees[i].tree.getNumberOfNodes() << " "
          //           << trees[j].tree.getNumberOfNodes()  << " ; "
          //           << time << "[µs]" << std::endl;
          dist = static_cast<double>(dist);
        }
        // std::cout << dist << " " << ttime << "[µs]" << std::endl;

        avg_time += ttime;
        avg_size += trees[r1].tree.getNumberOfNodes() + trees[r2].tree.getNumberOfNodes();
      }
      avg_size = avg_size / (nruns * 2);
      avg_time = avg_time / nruns;
      // std::cout << "------------\n";
      std::cout << "Avg size (random): " << avg_size << "\nAvg time (random): " << avg_time << std::endl;
      // for(auto r1: r1s){
      //   std::cout << r1 << ",";
      // }
      // std::cout << std::endl;
      // for(auto r2: r2s){
      //   std::cout << r2 << ",";
      // }
      // std::cout << std::endl;

      std::vector<std::vector<double>> dist_times(distanceMatrix.size(),std::vector<double>(distanceMatrix.size(),0));
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i) {
          // if(debugLevel_<3 and i % std::max(int(distanceMatrix.size() / 10), 1) == 0) {
          //   std::stringstream stream;
          //   stream << i << " / " << distanceMatrix.size();
          //   printMsg(stream.str());
          // }
          distanceMatrix[i][i] = 0.0;
          for(unsigned int j = i + 1; j < distanceMatrix[0].size(); ++j) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i,j) UNTIED() shared(distanceMatrix, dist_times, trees)
        {
#endif
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            // Execute
            if(baseModule_ == 0) {
              MergeTreeDistance mergeTreeDistance;
              mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
              mergeTreeDistance.setEpsilonTree1(epsilonTree1_);
              mergeTreeDistance.setEpsilonTree2(epsilonTree2_);
              mergeTreeDistance.setEpsilon2Tree1(epsilon2Tree1_);
              mergeTreeDistance.setEpsilon2Tree2(epsilon2Tree2_);
              mergeTreeDistance.setEpsilon3Tree1(epsilon3Tree1_);
              mergeTreeDistance.setEpsilon3Tree2(epsilon3Tree2_);
              mergeTreeDistance.setBranchDecomposition(branchDecomposition_);
              mergeTreeDistance.setParallelize(false);
              mergeTreeDistance.setPersistenceThreshold(persistenceThreshold_);
              mergeTreeDistance.setDebugLevel(std::min(debugLevel_, 2));
              mergeTreeDistance.setThreadNumber(this->threadNumber_);
              mergeTreeDistance.setNormalizedWasserstein(
                normalizedWasserstein_);
              mergeTreeDistance.setKeepSubtree(keepSubtree_);
              mergeTreeDistance.setDistanceSquaredRoot(distanceSquaredRoot_);
              mergeTreeDistance.setUseMinMaxPair(useMinMaxPair_);
              mergeTreeDistance.setPreprocess(false);
              // mergeTreeDistance.setSaveTree(true);
              mergeTreeDistance.setSaveTree(false);
              mergeTreeDistance.setCleanTree(true);
              mergeTreeDistance.setIsCalled(true);
              mergeTreeDistance.setPostprocess(false);
              mergeTreeDistance.setIsPersistenceDiagram(isPersistenceDiagram_);
              if(useDoubleInput_) {
                double const weight
                  = mixDistancesMinMaxPairWeight(isFirstInput);
                mergeTreeDistance.setMinMaxPairWeight(weight);
                mergeTreeDistance.setDistanceSquaredRoot(true);
              }
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> outputMatching;
              distanceMatrix[i][j] = mergeTreeDistance.execute<dataType>(
                trees[i], trees[j], outputMatching);
            } else if(baseModule_ == 1) {
              BranchMappingDistance branchDist;
              branchDist.setBaseMetric(branchMetric_);
              branchDist.setAssignmentSolver(assignmentSolverID_);
              branchDist.setSquared(distanceSquaredRoot_);
              branchDist.setEpsilonTree1(epsilonTree1_);
              branchDist.setEpsilonTree2(epsilonTree2_);
              branchDist.setEpsilon2Tree1(epsilon2Tree1_);
              branchDist.setEpsilon2Tree2(epsilon2Tree2_);
              branchDist.setEpsilon3Tree1(epsilon3Tree1_);
              branchDist.setEpsilon3Tree2(epsilon3Tree2_);
              branchDist.setPersistenceThreshold(persistenceThreshold_);
              branchDist.setPreprocess(false);
              // branchDist.setSaveTree(true);
              branchDist.setSaveTree(false);
              dataType dist = branchDist.execute<dataType>(trees[i], trees[j]);
              distanceMatrix[i][j] = static_cast<double>(dist);
            } else if(baseModule_ == 2) {
              PathMappingDistance pathDist;
              pathDist.setBaseMetric(pathMetric_);
              pathDist.setAssignmentSolver(assignmentSolverID_);
              pathDist.setSquared(distanceSquaredRoot_);
              pathDist.setComputeMapping(true);
              pathDist.setEpsilonTree1(epsilonTree1_);
              pathDist.setEpsilonTree2(epsilonTree2_);
              pathDist.setEpsilon2Tree1(epsilon2Tree1_);
              pathDist.setEpsilon2Tree2(epsilon2Tree2_);
              pathDist.setEpsilon3Tree1(epsilon3Tree1_);
              pathDist.setEpsilon3Tree2(epsilon3Tree2_);
              pathDist.setPersistenceThreshold(persistenceThreshold_);
              pathDist.setPreprocess(false);
              // pathDist.setSaveTree(true);
              pathDist.setSaveTree(false);
              pathDist.setlookahead(pathMappingLookahead_);
              if(vid_to_seg.size()==trees.size()) pathDist.setVidToSeg(vid_to_seg[i],vid_to_seg[j]);
              if(vid_to_seg.size()==trees.size()) pathDist.setValToSeg(val_to_seg[i],val_to_seg[j]);
              // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
              dataType dist = pathDist.execute<dataType>(trees[i], trees[j]);
              // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
              // auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
              // std::cout << i << " " << j << " ; "
              //           << trees[i].tree.getNumberOfNodes() << " "
              //           << trees[j].tree.getNumberOfNodes()  << " ; "
              //           << time << "[µs]" << std::endl;
              distanceMatrix[i][j] = static_cast<double>(dist);
            } else if(baseModule_ >= 3 && baseModule_ <= 5) {
              NaiveConstrMergeTreeEditDistance naiveEditDist;
              if(baseModule_==3) naiveEditDist.setconstraintType(0);
              if(baseModule_==4) naiveEditDist.setconstraintType(1);
              if(baseModule_==5) naiveEditDist.setconstraintType(2);
              naiveEditDist.setBaseMetric(pathMetric_);
              naiveEditDist.setAssignmentSolver(assignmentSolverID_);
              naiveEditDist.setSquared(distanceSquaredRoot_);
              naiveEditDist.setComputeMapping(true);
              naiveEditDist.setEpsilonTree1(epsilonTree1_);
              naiveEditDist.setEpsilonTree2(epsilonTree2_);
              naiveEditDist.setEpsilon2Tree1(epsilon2Tree1_);
              naiveEditDist.setEpsilon2Tree2(epsilon2Tree2_);
              naiveEditDist.setEpsilon3Tree1(epsilon3Tree1_);
              naiveEditDist.setEpsilon3Tree2(epsilon3Tree2_);
              naiveEditDist.setPersistenceThreshold(persistenceThreshold_);
              naiveEditDist.setPreprocess(false);
              // naiveEditDist.setSaveTree(true);
              naiveEditDist.setSaveTree(false);
              naiveEditDist.setlookahead(pathMappingLookahead_);
              if(vid_to_seg.size()==trees.size()) naiveEditDist.setVidToSeg(vid_to_seg[i],vid_to_seg[j]);
              if(vid_to_seg.size()==trees.size()) naiveEditDist.setValToSeg(val_to_seg[i],val_to_seg[j]);
              // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
              dataType dist = naiveEditDist.execute<dataType>(trees[i], trees[j]);
              // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
              // auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
              // std::cout << i << " " << j << " ; "
              //           << trees[i].tree.getNumberOfNodes() << " "
              //           << trees[j].tree.getNumberOfNodes()  << " ; "
              //           << time << "[µs]" << std::endl;
              distanceMatrix[i][j] = static_cast<double>(dist);
            }
            // distance matrix is symmetric
            distanceMatrix[j][i] = distanceMatrix[i][j];
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            auto single_dist_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
            // dist_times[i][j] = single_dist_time;
            if(debugLevel_>2) {
              std::stringstream stream;
              stream << i << "," << j << " / " << distanceMatrix.size() << "x" << distanceMatrix.size();
              printMsg(stream.str());
            }
#ifdef TTK_ENABLE_OPENMP
        } // end task
#endif
          } // end for j
          // if(debugLevel_>2) {
          //   std::stringstream stream;
          //   stream << i << " / " << distanceMatrix.size();
          //   printMsg(stream.str());
          // }
      } // end for i
      int num_dists = 0;
      int avg_matrix_time = 0;
      std::cout << "Number of threads: " << this->threadNumber_ << std::endl;
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i) {
        for(unsigned int j = i + 1; j < distanceMatrix[0].size(); ++j) {
          num_dists++;
        }
      }
      std::cout << "Number of instances: " << num_dists << std::endl;
      // for(unsigned int i = 0; i < distanceMatrix.size(); ++i) {
      //   for(unsigned int j = i + 1; j < distanceMatrix[0].size(); ++j) {
      //     avg_matrix_time += dist_times[i][j];
      //   }
      // }
      // std::cout << "Total dist time: " << avg_matrix_time/1000 << "ms" << std::endl;
      avg_matrix_time = avg_matrix_time/num_dists;
      // std::cout << "Average time: " << avg_matrix_time << "µs ; " << num_dists << " computations" << std::endl;
    }

  }; // MergeTreeDistanceMatrix class

} // namespace ttk
