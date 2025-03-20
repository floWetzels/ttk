/// \ingroup base
/// \class ttk::ComponentSize
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.03.2025
///
/// \brief TTK %componentSize processing package.
///
/// TODO

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>

#include <stack>

namespace ttk {
  class ComponentSize : virtual public Debug {

  public:
    struct Component {
      float center[3]{0, 0, 0};
      float size{0};
    };

    ComponentSize() {
      this->setDebugMsgPrefix("ComponentSize");
    };
    ~ComponentSize() override = default;

    int initializeComponents(
      std::unordered_map<int,Component>& components,
      const size_t n,
      const int* ids
     ) const {
      ttk::Timer t;
      this->printMsg("Initializing Components", 0, 0, ttk::debug::LineMode::REPLACE);

      std::vector<std::unordered_set<int>> unique_ids(this->threadNumber_);
      #pragma omp parallel num_threads(this->threadNumber_)
      {
        auto& ids_ = unique_ids[ omp_get_thread_num() ];

        #pragma omp for
        for (size_t i=0; i<n; i++) {
          ids_.insert(ids[i]);
        }
      }

      for(const auto& ids_ : unique_ids)
        for(const auto& id: ids_)
          components.emplace(id,Component());

      this->printMsg("Initializing Components (#"+std::to_string(components.size())+")", 1, t.getElapsedTime(), this->threadNumber_);
      return 1;
    };

    template <typename TT = ttk::AbstractTriangulation>
    int computeComponents(
      std::unordered_map<int,Component>& components,
      const size_t n,
      const int* ids,
      const TT *triangulation
    ){
      ttk::Timer t;
      this->printMsg("Computing Components", 0, 0, ttk::debug::LineMode::REPLACE);

      #pragma omp parallel num_threads(this->threadNumber_)
      {
        auto components_ = components;

        #pragma omp for
        for(size_t i=0; i<n; i++){
          const auto& id = ids[i];
          float coord[3];
          triangulation->getVertexPoint(i,coord[0],coord[1],coord[2]);
          auto& c = components_[id];
          c.size++;
          c.center[0]+=coord[0];
          c.center[1]+=coord[1];
          c.center[2]+=coord[2];
        }

        #pragma omp critical
        {
          for(const auto& it: components_){
            auto& c = components[it.first];
            c.size += it.second.size;
            c.center[0] += it.second.center[0];
            c.center[1] += it.second.center[1];
            c.center[2] += it.second.center[2];
          }
        }
      }

      for(auto& it: components){
        it.second.center[0] /= it.second.size;
        it.second.center[1] /= it.second.size;
        it.second.center[2] /= it.second.size;
      }

      this->printMsg("Computing Components", 1, t.getElapsedTime(), this->threadNumber_);
      return 1;
    };



  };
} // namespace ttk
