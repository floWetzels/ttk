/// \ingroup base
/// \class NaiveOneDegMergeTreeEditDistance
/// \author Florian Wetzels (wetzels@cs.uni-kl.de)
/// \date 2022.
///
/// This module defines the %NaiveOneDegMergeTreeEditDistance class that computes distances
/// between two merge trees.
///
/// \b Related \b publication \n
/// "A Deformation-based Edit Distance for Merge Trees" \n
/// Florian Wetzels, Christoph Garth. \n
/// TopoInVis 2022.

#pragma once

#include <set>
#include <vector>

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <vector>
#include <boost/dynamic_bitset.hpp>

// ttk common includes
#include "MergeTreeBase.h"
#include <AssignmentAuction.h>
#include <AssignmentExhaustive.h>
#include <AssignmentMunkres.h>
#include <Debug.h>
#include <FTMTree_MT.h>

namespace ttk {

  class NaiveOneDegMergeTreeEditDistance : virtual public Debug, public MergeTreeBase {

  private:
    int baseMetric_ = 0;
    int assignmentSolverID_ = 0;
    bool squared_ = false;
    bool computeMapping_ = false;
    unsigned int lookahead = 0;
    std::map<ftm::idNode,int> vid_to_seg1;
    std::map<ftm::idNode,int> vid_to_seg2;
    std::map<double,int> val_to_seg1;
    std::map<double,int> val_to_seg2;

    bool preprocess_ = true;
    bool saveTree_ = false;

    template <class dataType>
    inline dataType editCost_Persistence(int n1,
                                         int n2,
                                         ftm::FTMTree_MT *tree1,
                                         ftm::FTMTree_MT *tree2) {
      if(!vid_to_seg1.empty() and !vid_to_seg2.empty()
         and n1>0 and tree1->getNumberOfChildren(n1)==0
         and n2>0 and tree2->getNumberOfChildren(n2)==0){
        auto val1 = tree1->getValue<dataType>(n1);
        auto val2 = tree2->getValue<dataType>(n2);
        auto s1 = val_to_seg1[val1];
        auto s2 = val_to_seg2[val2];
        if(s1!=s2){
          dataType b1 = tree1->getValue<dataType>(n1);
          dataType d1 = tree1->getValue<dataType>(tree1->getParent(n1));
          dataType b2 = tree2->getValue<dataType>(n2);
          dataType d2 = tree2->getValue<dataType>(tree2->getParent(n2));
          dataType del1 = (d1 > b1) ? (d1 - b1) : (b1 - d1);
          dataType del2 = (d2 > b2) ? (d2 - b2) : (b2 - d2);
          return (del1+del2)*2;
        }
      }
      dataType d;
      if(n1 < 0) {
        dataType b1 = tree2->getValue<dataType>(n2);
        dataType d1 = tree2->getValue<dataType>(tree2->getParent(n2));
        d = (d1 > b1) ? (d1 - b1) : (b1 - d1);
      } else if(n2 < 0) {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(tree1->getParent(n1));
        d = (d1 > b1) ? (d1 - b1) : (b1 - d1);
      } else {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(tree1->getParent(n1));
        dataType b2 = tree2->getValue<dataType>(n2);
        dataType d2 = tree2->getValue<dataType>(tree2->getParent(n2));
        dataType dist1 = (d1 > b1) ? (d1 - b1) : (b1 - d1);
        dataType dist2 = (d2 > b2) ? (d2 - b2) : (b2 - d2);
        d = (dist1 > dist2) ? (dist1 - dist2) : (dist2 - dist1);
      }
      return squared_ ? d * d : d;
    }

    template <class dataType>
    void traceMapping_path(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      int curr1,
      int curr2,
      std::vector<std::vector<int>> &predecessors1,
      std::vector<std::vector<int>> &predecessors2,
      int depth1,
      int depth2,
      std::vector<dataType> &memT,
      std::vector<std::pair<dataType,std::vector<std::pair<ftm::idNode, ftm::idNode>>>> &memLA,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &mapping) {

      //===============================================================================
      // If both trees not empty, find optimal edit operation
      std::vector<ftm::idNode> children1;
      tree1->getChildren(curr1, children1);
      std::vector<ftm::idNode> children2;
      tree2->getChildren(curr2, children2);
      // int parent1 = predecessors1[curr1][predecessors1[curr1].size() - l1];
      // int parent2 = predecessors2[curr2][predecessors2[curr2].size() - l2];

      size_t const nn1 = tree1->getNumberOfNodes();
      size_t const nn2 = tree2->getNumberOfNodes();
      size_t const dim1 = 1;
      size_t const dim2 = (nn1 + 1) * dim1;
      size_t const dim3 = (depth1 + 1) * dim2;
      size_t const dim4 = (nn2 + 1) * dim3;

      //---------------------------------------------------------------------------
      // If both trees only have one branch, return edit cost between
      // the two branches
      if(tree1->getNumberOfChildren(curr1) == 0
         and tree2->getNumberOfChildren(curr2) == 0) {
        mapping.emplace_back(
          curr1,curr2,editCost_Persistence<dataType>(curr1,curr2,tree1,tree2));
        return;
      }
      //---------------------------------------------------------------------------
      // If both trees have more than one branch, try all decompositions
      // of both trees
      else {
        //-----------------------------------------------------------------------
        // Try all possible main branches of first tree (child1_mb) and
        // all possible main branches of second tree (child2_mb) Then
        // try all possible matchings of subtrees
        if(tree1->getNumberOfChildren(curr1) <= 2
           && tree2->getNumberOfChildren(curr2) <= 2) {
          int const child11 = children1.size()>0 ? children1[0] : nn1;
          int const child12 = children1.size()>1 ? children1[1] : nn1;
          int const child21 = children2.size()>0 ? children2[0] : nn2;
          int const child22 = children2.size()>1 ? children2[1] : nn2;
          if(memT[curr1 + curr2 * (nn1+1)]
             == memT[child11 + child21 * (nn1+1)]
                  + memT[child12 + child22 * (nn1+1)]
                  + editCost_Persistence<dataType>(
                    curr1, curr2, tree1, tree2)) {
            mapping.emplace_back(
              curr1,curr2,editCost_Persistence<dataType>(curr1,curr2,tree1,tree2));
            if(child11<nn1 && child21<nn2){
              traceMapping_path(tree1, tree2, child11, child21,
                                predecessors1, predecessors2, depth1, depth2,
                                memT, memLA, mapping);
            }
            if(child12<nn1 && child22<nn2){
              traceMapping_path(tree1, tree2, child12, child22,
                                predecessors1, predecessors2, depth1, depth2,
                                memT, memLA, mapping);
            }
            return;
          }
          if(memT[curr1 + curr2 * (nn1+1)]
             == memT[child11 + child22 * (nn1+1)]
                  + memT[child12 + child21 * (nn1+1)]
                  + editCost_Persistence<dataType>(
                    curr1, curr2, tree1, tree2)) {
            mapping.emplace_back(
              curr1,curr2,editCost_Persistence<dataType>(curr1,curr2,tree1,tree2));
            if(child11<nn1 && child22<nn2){
              traceMapping_path(tree1, tree2, child11, child22,
                                predecessors1, predecessors2, depth1, depth2,
                                memT, memLA, mapping);
            }
            if(child12<nn1 && child21<nn2){
              traceMapping_path(tree1, tree2, child12, child21,
                                predecessors1, predecessors2, depth1, depth2,
                                memT, memLA, mapping);
            }
            return;
          }
        } else {
          auto f = [&](int r, int c) {
            size_t const c1
              = r < tree1->getNumberOfChildren(curr1) ? children1[r] : nn1;
            size_t const c2
              = c < tree2->getNumberOfChildren(curr2) ? children2[c] : nn2;
            int const l1_ = c1 == nn1 ? 0 : 1;
            int const l2_ = c2 == nn2 ? 0 : 1;
            return memT[c1 + c2 * (nn1+1)];
          };
          int size = std::max(tree1->getNumberOfChildren(curr1),
                              tree2->getNumberOfChildren(curr2))
                     + 1;
          auto costMatrix = std::vector<std::vector<dataType>>(
            size, std::vector<dataType>(size, 0));
          std::vector<MatchingType> matching;
          for(int r = 0; r < size; r++) {
            for(int c = 0; c < size; c++) {
              costMatrix[r][c] = f(r, c);
            }
          }

          AssignmentSolver<dataType> *assignmentSolver;
          AssignmentExhaustive<dataType> solverExhaustive;
          AssignmentMunkres<dataType> solverMunkres;
          AssignmentAuction<dataType> solverAuction;
          switch(assignmentSolverID_) {
            case 1:
              solverExhaustive = AssignmentExhaustive<dataType>();
              assignmentSolver = &solverExhaustive;
              break;
            case 2:
              solverMunkres = AssignmentMunkres<dataType>();
              assignmentSolver = &solverMunkres;
              break;
            case 0:
            default:
              solverAuction = AssignmentAuction<dataType>();
              assignmentSolver = &solverAuction;
          }
          assignmentSolver->setInput(costMatrix);
          assignmentSolver->setBalanced(true);
          assignmentSolver->run(matching);
          dataType d_ = editCost_Persistence<dataType>(
            curr1, curr2, tree1, tree2);
          for(auto m : matching)
            d_ += std::get<2>(m);
          if(memT[curr1 + curr2 * (nn1+1)] == d_) {
            mapping.emplace_back(
              curr1,curr2,editCost_Persistence<dataType>(curr1,curr2,tree1,tree2));
            for(auto m : matching) {
              int n1 = std::get<0>(m) < tree1->getNumberOfChildren(curr1)
                         ? children1[std::get<0>(m)]
                         : -1;
              int n2 = std::get<1>(m) < tree2->getNumberOfChildren(curr2)
                         ? children2[std::get<1>(m)]
                         : -1;
              if(n1 >= 0 && n2 >= 0)
                traceMapping_path(tree1, tree2, n1, n2, predecessors1,
                                  predecessors2, depth1, depth2, memT, memLA, mapping);
            }
            return;
          }
        }    
        //-----------------------------------------------------------------------
        // Try to look-ahead
        if(memLA.size()>0 and memLA[curr1+curr2*nn1].first>=0){
          dataType d_ = editCost_Persistence<dataType>(
              curr1, curr2, tree1, tree2) 
              +memLA[curr1+curr2*nn1].first;
          if(memT[curr1 + curr2 * (nn1+1)] == d_){
            mapping.emplace_back(
                curr1,curr2,editCost_Persistence<dataType>(curr1,curr2,tree1,tree2));
            for(auto m : memLA[curr1+curr2*nn1].second){
              traceMapping_path(tree1, tree2, m.first, m.second, predecessors1,
                                    predecessors2, depth1, depth2, memT, memLA, mapping);
            }
            return;
          }
        }
      }
    }

  public:
    NaiveOneDegMergeTreeEditDistance() {
      this->setDebugMsgPrefix(
        "MergeTreeDistance"); // inherited from Debug: prefix will be printed at
                              // the beginning of every msg
    }
    ~NaiveOneDegMergeTreeEditDistance() override = default;

    void setBaseMetric(int m) {
      baseMetric_ = m;
    }

    void setAssignmentSolver(int assignmentSolver) {
      assignmentSolverID_ = assignmentSolver;
    }

    void setSquared(bool s) {
      squared_ = s;
    }

    void setComputeMapping(bool m) {
      computeMapping_ = m;
    }

    void setlookahead(int l) {
      lookahead = l;
    }

    void setPreprocess(bool p) {
      preprocess_ = p;
    }

    void setSaveTree(bool save) {
      saveTree_ = save;
    }

    void setVidToSeg(std::map<ftm::idNode,int> m1, std::map<ftm::idNode,int> m2) {
      vid_to_seg1 = m1;
      vid_to_seg2 = m2;
    }

    void setValToSeg(std::map<double,int> m1, std::map<double,int> m2) {
      val_to_seg1 = m1;
      val_to_seg2 = m2;
    }

    template <class dataType>
    dataType computeDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode,double>> *outputMatching) {

      // compute preorder of both trees (necessary for bottom-up dynamic
      // programming)

      std::vector<std::vector<int>> predecessors1(tree1->getNumberOfNodes());
      std::vector<std::vector<int>> predecessors2(tree2->getNumberOfNodes());
      int const rootID1 = tree1->getRoot();
      int const rootID2 = tree2->getRoot();
      std::vector<int> preorder1(tree1->getNumberOfNodes());
      std::vector<int> preorder2(tree2->getNumberOfNodes());
      std::vector<dataType> scalarDepth1(tree1->getNumberOfNodes());
      std::vector<dataType> totalPersistence1(tree1->getNumberOfNodes());
      std::vector<dataType> scalarDepth2(tree2->getNumberOfNodes());
      std::vector<dataType> totalPersistence2(tree2->getNumberOfNodes());

      int depth1 = 0;
      int depth2 = 0;
      std::stack<int> stack;
      std::stack<int> postorderstack;
      stack.push(rootID1);
      int count = tree1->getNumberOfNodes() - 1;
      while(!stack.empty()) {
        int const nIdx = stack.top();
        stack.pop();
        postorderstack.push(nIdx);
        preorder1[count] = nIdx;
        count--;
        depth1 = std::max((int)predecessors1[nIdx].size(), depth1);
        std::vector<ftm::idNode> children;
        tree1->getChildren(nIdx, children);
        for(int const cIdx : children) {
          stack.push(cIdx);
          predecessors1[cIdx].reserve(predecessors1[nIdx].size() + 1);
          predecessors1[cIdx].insert(predecessors1[cIdx].end(),
                                     predecessors1[nIdx].begin(),
                                     predecessors1[nIdx].end());
          predecessors1[cIdx].push_back(nIdx);
        }
      }
      while(!postorderstack.empty()) {
        int const nIdx = postorderstack.top();
        postorderstack.pop();
        dataType total = 0;
        dataType longest = 0;
        std::vector<ftm::idNode> children;
        tree1->getChildren(nIdx, children);
        auto nv = tree1->getValue<dataType>(nIdx);
        for(int const cIdx : children) {
          auto cv = tree1->getValue<dataType>(cIdx);
          auto pers = nv > cv ? nv-cv : cv-nv;
          total += pers + totalPersistence1[cIdx];
          longest = std::max(longest,pers + scalarDepth1[cIdx]);
        }
        totalPersistence1[nIdx] = total;
        scalarDepth1[nIdx] = longest;
      }
      postorderstack = std::stack<int>();
      stack.push(rootID2);
      count = tree2->getNumberOfNodes() - 1;
      while(!stack.empty()) {
        int const nIdx = stack.top();
        stack.pop();
        postorderstack.push(nIdx);
        preorder2[count] = nIdx;
        count--;
        depth2 = std::max((int)predecessors2[nIdx].size(), depth2);
        std::vector<ftm::idNode> children;
        tree2->getChildren(nIdx, children);
        for(int const cIdx : children) {
          stack.push(cIdx);
          predecessors2[cIdx].reserve(predecessors2[nIdx].size() + 1);
          predecessors2[cIdx].insert(predecessors2[cIdx].end(),
                                     predecessors2[nIdx].begin(),
                                     predecessors2[nIdx].end());
          predecessors2[cIdx].push_back(nIdx);
        }
      }
      while(!postorderstack.empty()) {
        int const nIdx = postorderstack.top();
        postorderstack.pop();
        dataType total = 0;
        dataType longest = 0;
        std::vector<ftm::idNode> children;
        tree2->getChildren(nIdx, children);
        auto nv = tree2->getValue<dataType>(nIdx);
        for(int const cIdx : children) {
          auto cv = tree2->getValue<dataType>(cIdx);
          auto pers = nv > cv ? nv-cv : cv-nv;
          total += pers + totalPersistence2[cIdx];
          longest = std::max(longest,pers + scalarDepth2[cIdx]);
        }
        totalPersistence2[nIdx] = total;
        scalarDepth2[nIdx] = longest;
      }

      dataType tP1 = totalPersistence1[rootID1];
      dataType tP2 = totalPersistence2[rootID2];
      dataType lP1 = scalarDepth1[rootID1];
      dataType lP2 = scalarDepth2[rootID2];
      // dataType globalLowerBound = tP1>tP2 ? tP1-tP2 : tP2-tP1;
      dataType globalUpperBound = lP1>lP2 ? lP1-lP2 : lP2-lP1;
      globalUpperBound += tP1-lP1 + tP2-lP2;

      if(lookahead>0){
        auto lookahead_tmp = lookahead;
        auto computeMapping_tmp = computeMapping_;
        auto preprocess_tmp = preprocess_;
        lookahead = 0;
        computeMapping_ = false;

        auto dist = computeDistance<dataType>(tree1,tree2,outputMatching);
        globalUpperBound = dist;

        preprocess_ = preprocess_tmp;
        computeMapping_ = computeMapping_tmp;
        lookahead = lookahead_tmp;
      }

      // initialize memoization tables

      size_t nn1 = tree1->getNumberOfNodes();
      size_t nn2 = tree2->getNumberOfNodes();

      std::vector<dataType> memT((nn1 + 1) * (nn2 + 1));
      std::vector<std::pair<dataType,std::vector<std::pair<ftm::idNode, ftm::idNode>>>> memLA;
      if(lookahead>0 and computeMapping_){
        memLA.resize((nn1) * (nn2),std::make_pair(-1.0,std::vector<std::pair<ftm::idNode, ftm::idNode>>()));
      }

      //---------------------------------------------------------------------------
      // fill memoization table with base cases (subtree deletions)
      memT[nn1 + nn2 * (nn1+1)] = 0;
      for(size_t i = 0; i < nn1; i++) {
        if(i==rootID1) continue;
        int curr1 = preorder1[i];
        std::vector<ftm::idNode> children1;
        tree1->getChildren(curr1, children1);

        //-----------------------------------------------------------------------
        // Delete curr path and full subtree rooted in path
        memT[curr1 + nn2 * (nn1+1)]
          = editCost_Persistence<dataType>(
            curr1, -1, tree1, tree2);
        for(auto child1 : children1) {
          memT[curr1 + nn2 * (nn1+1)]
            += memT[child1 + nn2 * (nn1+1)];
        }
      }
      for(size_t j = 0; j < nn2; j++) {
        if(j==rootID2) continue;
        int curr2 = preorder2[j];
        std::vector<ftm::idNode> children2;
        tree2->getChildren(curr2, children2);

        //-----------------------------------------------------------------------
        // Delete curr path and full subtree rooted in path
        memT[nn1 + curr2 * (nn1+1)]
          = editCost_Persistence<dataType>(
            -1, curr2, tree1, tree2);
        for(auto child2 : children2) {
          memT[nn1 + curr2 * (nn1+1)]
            += memT[nn1 + child2 * (nn1+1)];
        }
      }

      //---------------------------------------------------------------------------
      // fill memoization table with cases for two non-empty subtrees
      for(size_t i = 0; i < nn1; i++) {
        if(i==rootID1) continue;
        int curr1 = preorder1[i];
        std::vector<ftm::idNode> children1;
        tree1->getChildren(curr1, children1);
        for(size_t j = 0; j < nn2; j++) {
          if(j==rootID2) continue;
          int curr2 = preorder2[j];
          std::vector<ftm::idNode> children2;
          tree2->getChildren(curr2, children2);

          //===============================================================================
          // check if lookahead is necessary
          bool hasSaddleChildren1 = false;
          bool hasSaddleChildren2 = false;
          for (auto c : children1){
            if(!tree1->isLeaf(c)){
              hasSaddleChildren1 = true;
              break;
            }
          }
          for (auto c : children2){
            if(!tree2->isLeaf(c)){
              hasSaddleChildren2 = true;
              break;
            }
          }
          tP1 = totalPersistence1[curr1];
          tP2 = totalPersistence2[curr2];
          dataType localLowerBound = tP1>tP2 ? tP1-tP2 : tP2-tP1;
          bool useLookahead = lookahead>0
                              and !children1.empty() and !children2.empty()
                              and hasSaddleChildren1 and hasSaddleChildren2
                              and localLowerBound<globalUpperBound;

          //===============================================================================
          // compute optimal lookahead mapping between children if necessary
          if(useLookahead){  
            std::stack<std::pair<ftm::idNode,ftm::idNode>> s;
            std::unordered_map<ftm::idNode,ftm::idNode> right;
            std::unordered_map<ftm::idNode,ftm::idNode> down;
            for (ftm::idNode ci=0; ci<children1.size()-1; ci++){
              s.push(std::make_pair(children1[ci],children1[ci+1]));
            }
            s.push(std::make_pair(children1.back(),children1.back()));
            while(not s.empty()){
              auto cn = s.top().first;
              auto cr = s.top().second;
              s.pop();
              if(predecessors1[cn].size()-predecessors1[curr1].size()>lookahead+1){
                continue;
              }
              right[cn] = cr;
              if(tree1->getNumberOfChildren(cn)>0){
                std::vector<unsigned int> children_cn;
                tree1->getChildren(cn, children_cn);
                down[cn] = children_cn.front();
                for(ftm::idNode ci=0; ci<children_cn.size()-1; ci++){
                  s.push(std::make_pair(children_cn[ci],children_cn[ci+1]));
                }
                s.push(std::make_pair(children_cn.back(),cn==cr?children_cn.back():cr));
              }
            }
            std::vector<std::pair<dataType,std::vector<unsigned int>>> cases1;
            std::vector<std::tuple<ftm::idNode,dataType,std::vector<unsigned int>>> worklist;
            cases1.reserve(right.size());
            worklist.reserve(right.size());
            worklist.push_back(std::make_tuple(children1[0],0,std::vector<unsigned int>()));
            while (!worklist.empty()){
              auto curr_tuple = worklist.back();
              auto cn = std::get<0>(curr_tuple);
              auto deleted_cost = std::get<1>(curr_tuple);
              auto kept_nodes = std::get<2>(curr_tuple);
              worklist.pop_back();
              auto p = predecessors1[cn].back();
              auto l = predecessors1[cn].size()-predecessors1[curr1].size();
              if (l<=lookahead and tree1->getNumberOfChildren(cn)>0){
                auto deletion_cost = editCost_Persistence<dataType>(cn,-1,tree1,tree2);
                worklist.push_back(std::make_tuple(down[cn],deleted_cost+deletion_cost,kept_nodes));
              }
              std::vector<unsigned int> kept_nodes_ = kept_nodes;
              kept_nodes_.push_back(cn);
              if(right[cn]==cn){
                if(kept_nodes_.size()>1 && deleted_cost>0) cases1.push_back(std::make_pair(deleted_cost,kept_nodes_));
              }
              else{
                worklist.push_back(std::make_tuple(right[cn],deleted_cost,kept_nodes_));
              }
            }

            s = std::stack<std::pair<ftm::idNode,ftm::idNode>>();
            right = std::unordered_map<ftm::idNode,ftm::idNode>();
            down = std::unordered_map<ftm::idNode,ftm::idNode>();
            for (ftm::idNode ci=0; ci<children2.size()-1; ci++){
              s.push(std::make_pair(children2[ci],children2[ci+1]));
            }
            s.push(std::make_pair(children2.back(),children2.back()));
            while(not s.empty()){
              ftm::idNode cn = s.top().first;
              ftm::idNode cr = s.top().second;
              s.pop();
              if(predecessors2[cn].size()-predecessors2[curr2].size()>lookahead+1){
                continue;
              }
              right[cn] = cr;
              if(tree2->getNumberOfChildren(cn)>0){
                std::vector<unsigned int> children_cn;
                tree2->getChildren(cn, children_cn);
                down[cn] = children_cn.front();
                for(ftm::idNode ci=0; ci<children_cn.size()-1; ci++){
                  s.push(std::make_pair(children_cn[ci],children_cn[ci+1]));
                }
                s.push(std::make_pair(children_cn.back(),cn==cr?children_cn.back():cr));
              }
            }
            std::vector<std::pair<dataType,std::vector<unsigned int>>> cases2;
            worklist = std::vector<std::tuple<ftm::idNode,dataType,std::vector<unsigned int>>>();
            cases2.reserve(right.size());
            worklist.reserve(right.size());
            worklist.push_back(std::make_tuple(children2[0],0,std::vector<unsigned int>()));
            while (!worklist.empty()){
              auto curr_tuple = worklist.back();
              auto cn = std::get<0>(curr_tuple);
              auto deleted_cost = std::get<1>(curr_tuple);
              auto kept_nodes = std::get<2>(curr_tuple);
              worklist.pop_back();
              auto p = predecessors2[cn].back();
              auto l = predecessors2[cn].size()-predecessors2[curr2].size();
              if (l<=lookahead and tree2->getNumberOfChildren(cn)>0){
                auto deletion_cost = editCost_Persistence<dataType>(-1,cn,tree1,tree2);
                worklist.push_back(std::make_tuple(down[cn],deleted_cost+deletion_cost,kept_nodes));
              }
              std::vector<unsigned int> kept_nodes_ = kept_nodes;
              kept_nodes_.push_back(cn);
              if(right[cn]==cn){
                if(kept_nodes_.size()>1 && deleted_cost>0) cases2.push_back(std::make_pair(deleted_cost,kept_nodes_));
              }
              else{
                worklist.push_back(std::make_tuple(right[cn],deleted_cost,kept_nodes_));
              }
            }

            dataType opt_case_cost = std::numeric_limits<dataType>::max();
            std::vector<std::pair<ftm::idNode,ftm::idNode>> opt_case;
            for (auto case1 : cases1){
              auto delete_costs1 = std::get<0>(case1);
              if(delete_costs1>globalUpperBound) continue;
              auto actual_children1 = std::get<1>(case1);
              dataType tP1_ = tP1-delete_costs1;
              dataType tP1__ = totalPersistence1[rootID1]-tP1;
              for (auto case2 : cases2){
                auto delete_costs2 = std::get<0>(case2);
                if(delete_costs1+delete_costs2>globalUpperBound) continue;
                dataType case_cost = delete_costs1+delete_costs2;
                dataType tP2_ = tP2-delete_costs2;
                dataType tP2__ = totalPersistence2[rootID2]-tP2;
                // dataType localBound = (tP1_>tP2_ ? tP1_-tP2_ : tP2_-tP1_) + case_cost;
                dataType localBound = (tP1__>tP2__ ? tP1__-tP2__ : tP2__-tP1__)
                                    + (tP1_>tP2_ ? tP1_-tP2_ : tP2_-tP1_) 
                                    + case_cost;
                if(localBound>globalUpperBound) continue;
                auto actual_children2 = std::get<1>(case2);

                auto f = [&](unsigned int r, unsigned int c) {
                  size_t const c1 = r < actual_children1.size()
                                      ? actual_children1[r]
                                      : nn1;
                  size_t const c2 = c < actual_children2.size()
                                      ? actual_children2[c]
                                      : nn2;
                  int const l1_ = c1 == nn1 ? 0 : 1;
                  int const l2_ = c2 == nn2 ? 0 : 1;
                  return memT[c1 + c2 * (nn1+1)];
                };
                int size = std::max(actual_children1.size(),actual_children2.size()) + 1;
                auto costMatrix = std::vector<std::vector<dataType>>(
                  size, std::vector<dataType>(size, 0));
                std::vector<MatchingType> matching;
                for(int r = 0; r < size; r++) {
                  for(int c = 0; c < size; c++) {
                    costMatrix[r][c] = f(r, c);
                  }
                }

                AssignmentSolver<dataType> *assignmentSolver;
                AssignmentExhaustive<dataType> solverExhaustive;
                AssignmentMunkres<dataType> solverMunkres;
                AssignmentAuction<dataType> solverAuction;
                switch(assignmentSolverID_) {
                  case 1:
                    solverExhaustive = AssignmentExhaustive<dataType>();
                    assignmentSolver = &solverExhaustive;
                    break;
                  case 2:
                    solverMunkres = AssignmentMunkres<dataType>();
                    assignmentSolver = &solverMunkres;
                    break;
                  case 0:
                  default:
                    solverAuction = AssignmentAuction<dataType>();
                    assignmentSolver = &solverAuction;
                }
                assignmentSolver->setInput(costMatrix);
                assignmentSolver->setBalanced(true);
                assignmentSolver->run(matching);
                for(auto m : matching)
                  case_cost += std::get<2>(m);
                opt_case_cost = std::min(opt_case_cost, case_cost);
                if(opt_case_cost==case_cost){
                  opt_case = std::vector<std::pair<ftm::idNode,ftm::idNode>>();
                  for(auto m : matching){
                    ftm::idNode m1 = std::get<0>(m);
                    ftm::idNode m2 = std::get<1>(m);
                    if(m1<actual_children1.size() and m2<actual_children2.size()){
                      opt_case.push_back(std::make_pair(actual_children1[m1],actual_children2[m2]));
                    }
                  }
                }
              }
            }

            // memT[curr1 + curr2 * (nn1+1)] = opt_case_cost;
            memLA[curr1+curr2*nn1] = std::make_pair(opt_case_cost,opt_case);
          }

          //===============================================================================
          // normal recursions for two non-empty subtrees

          //---------------------------------------------------------------------------
          // If both trees only have one edge, return edit cost between
          // the two edges
          if(tree1->getNumberOfChildren(curr1) == 0
              and tree2->getNumberOfChildren(curr2) == 0) {
            memT[curr1 + curr2 * (nn1+1)]
              = editCost_Persistence<dataType>(
                curr1, curr2, tree1, tree2);
          }
          //---------------------------------------------------------------------------
          // If both trees have more than one edge, try matching cases and deletion cases
          else {
            dataType d = std::numeric_limits<dataType>::max();
            //-----------------------------------------------------------------------
            // Try matching both root edges and then find optimal assignment between subtrees
            if(tree1->getNumberOfChildren(curr1) <= 2
                && tree2->getNumberOfChildren(curr2) <= 2) {
              int const child11 = children1.size()>0 ? children1[0] : nn1;
              int const child12 = children1.size()>1 ? children1[1] : nn1;
              int const child21 = children2.size()>0 ? children2[0] : nn2;
              int const child22 = children2.size()>1 ? children2[1] : nn2;
              d = std::min<dataType>(
                d, memT[child11 + child21 * (nn1+1)]
                      + memT[child12 + child22 * (nn1+1)]
                      + editCost_Persistence<dataType>(
                        curr1, curr2, tree1, tree2));
              d = std::min<dataType>(
                d, memT[child11 + child22 * (nn1+1)]
                      + memT[child12 + child21 * (nn1+1)]
                      + editCost_Persistence<dataType>(
                        curr1, curr2, tree1, tree2));
            } else {
              auto f = [&](int r, int c) {
                size_t const c1 = r < tree1->getNumberOfChildren(curr1)
                                    ? children1[r]
                                    : nn1;
                size_t const c2 = c < tree2->getNumberOfChildren(curr2)
                                    ? children2[c]
                                    : nn2;
                int const l1_ = c1 == nn1 ? 0 : 1;
                int const l2_ = c2 == nn2 ? 0 : 1;
                return memT[c1 + c2 * (nn1+1)];
              };
              int size = std::max(tree1->getNumberOfChildren(curr1),
                                  tree2->getNumberOfChildren(curr2))
                          + 1;
              auto costMatrix = std::vector<std::vector<dataType>>(
                size, std::vector<dataType>(size, 0));
              std::vector<MatchingType> matching;
              for(int r = 0; r < size; r++) {
                for(int c = 0; c < size; c++) {
                  costMatrix[r][c] = f(r, c);
                }
              }

              AssignmentSolver<dataType> *assignmentSolver;
              AssignmentExhaustive<dataType> solverExhaustive;
              AssignmentMunkres<dataType> solverMunkres;
              AssignmentAuction<dataType> solverAuction;
              switch(assignmentSolverID_) {
                case 1:
                  solverExhaustive = AssignmentExhaustive<dataType>();
                  assignmentSolver = &solverExhaustive;
                  break;
                case 2:
                  solverMunkres = AssignmentMunkres<dataType>();
                  assignmentSolver = &solverMunkres;
                  break;
                case 0:
                default:
                  solverAuction = AssignmentAuction<dataType>();
                  assignmentSolver = &solverAuction;
              }
              assignmentSolver->setInput(costMatrix);
              assignmentSolver->setBalanced(true);
              assignmentSolver->run(matching);
              dataType d_ = editCost_Persistence<dataType>(
                curr1, curr2, tree1, tree2);
              for(auto m : matching)
                d_ += std::get<2>(m);
              d = std::min(d, d_);
            }
            
            //-----------------------------------------------------------------------
            // Try to match both root edges and optimal look-ahead assignment between subtrees
            if(useLookahead){
              dataType case_cost = editCost_Persistence<dataType>(curr1,curr2,tree1,tree2);
              case_cost += memLA[curr1 + curr2 * nn1].first;
              d = std::min(d, case_cost);
            }
            memT[curr1 + curr2 * (nn1+1)] = d;
          }
        }
      }

      std::vector<ftm::idNode> children1;
      tree1->getChildren(rootID1, children1);
      std::vector<ftm::idNode> children2;
      tree2->getChildren(rootID2, children2);

      dataType res
        = memT[children1[0] + children2[0] * (nn1+1)];

      if(computeMapping_ && outputMatching) {

        outputMatching->clear();
        traceMapping_path(tree1, tree2, children1[0], children2[0],
                          predecessors1, predecessors2, depth1, depth2, memT, memLA,
                          *outputMatching);

      }

      return squared_ ? std::sqrt(res) : res;
    }

    template <class dataType>
    dataType execute(ftm::MergeTree<dataType> &mTree1,
                     ftm::MergeTree<dataType> &mTree2,
                     std::vector<std::tuple<ftm::idNode, ftm::idNode,double>> *outputMatching) {

      ftm::MergeTree<dataType> mTree1Copy;
      ftm::MergeTree<dataType> mTree2Copy;
      if(saveTree_) {
        mTree1Copy = ftm::copyMergeTree<dataType>(mTree1);
        mTree2Copy = ftm::copyMergeTree<dataType>(mTree2);
      }
      ftm::MergeTree<dataType> &mTree1Int = (saveTree_ ? mTree1Copy : mTree1);
      ftm::MergeTree<dataType> &mTree2Int = (saveTree_ ? mTree2Copy : mTree2);
      ftm::FTMTree_MT *tree1 = &(mTree1Int.tree);
      ftm::FTMTree_MT *tree2 = &(mTree2Int.tree);

      // optional preprocessing
      if(preprocess_) {
        treesNodeCorr_.resize(2);
        preprocessingPipeline<dataType>(
          mTree1Int, epsilonTree1_, epsilon2Tree1_, epsilon3Tree1_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[0],
          true, true);
        preprocessingPipeline<dataType>(
          mTree2Int, epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[1],
          true, true);
      }

      tree1 = &(mTree1Int.tree);
      tree2 = &(mTree2Int.tree);

      return computeDistance<dataType>(tree1, tree2, outputMatching);
    }

    template <class dataType>
    dataType computeDistance(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2) {
      return computeDistance<dataType>(
        tree1, tree2,
        (std::vector<std::tuple<ftm::idNode, ftm::idNode,double>> *)nullptr);
    }

    template <class dataType>
    dataType execute(ftm::MergeTree<dataType> &mTree1,
                     ftm::MergeTree<dataType> &mTree2) {

      ftm::MergeTree<dataType> mTree1Copy;
      ftm::MergeTree<dataType> mTree2Copy;
      if(saveTree_) {
        mTree1Copy = ftm::copyMergeTree<dataType>(mTree1);
        mTree2Copy = ftm::copyMergeTree<dataType>(mTree2);
      }
      ftm::MergeTree<dataType> &mTree1Int = (saveTree_ ? mTree1Copy : mTree1);
      ftm::MergeTree<dataType> &mTree2Int = (saveTree_ ? mTree2Copy : mTree2);
      ftm::FTMTree_MT *tree1 = &(mTree1Int.tree);
      ftm::FTMTree_MT *tree2 = &(mTree2Int.tree);

      // optional preprocessing
      if(preprocess_) {
        treesNodeCorr_.resize(2);
        preprocessingPipeline<dataType>(
          mTree1Int, epsilonTree1_, epsilon2Tree1_, epsilon3Tree1_, false,
          useMinMaxPair_, cleanTree_, treesNodeCorr_[0], true, true);
        preprocessingPipeline<dataType>(
          mTree2Int, epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_, false,
          useMinMaxPair_, cleanTree_, treesNodeCorr_[1], true, true);
      }

      tree1 = &(mTree1Int.tree);
      tree2 = &(mTree2Int.tree);

      return computeDistance<dataType>(
        tree1, tree2,
        (std::vector<std::tuple<ftm::idNode, ftm::idNode,double>> *)nullptr);
    }
  };
} // namespace ttk
