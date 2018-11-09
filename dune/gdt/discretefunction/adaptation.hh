// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_ADAPTATION_HH
#define DUNE_GDT_DISCRETEFUNCTION_ADAPTATION_HH

#include <list>

#include <dune/common/dynvector.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/conversion.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * \todo Only depend on the grid type, use type erasure to allow to add arbitrary discrete functions, possibly implement
 *       this by element functors which are used in a grid walker.
 */
template <class V, class GV, size_t r = 1, size_t rC = 1, class RF = double>
class AdaptationHelper
{
  using ThisType = AdaptationHelper<V, GV, r, rC, RF>;

public:
  using DiscreteFunctionType = DiscreteFunction<V, GV, r, rC, RF>;
  using SpaceType = SpaceInterface<GV, r, rC, RF>;
  using G = typename GV::Grid;

  AdaptationHelper(G& grd)
    : grid_(grd)
    , data_()
  {
  }

  ThisType& append(SpaceType& space, DiscreteFunctionType& discrete_function)
  {
    data_.emplace_back(space,
                       discrete_function,
                       PersistentContainer<G, DynamicVector<RF>>(grid_, 0),
                       discrete_function.local_discrete_function());
    return *this;
  }

  void pre_adapt(const bool pre_adapt_grid = true)
  {
    auto grid_view = grid_.leafGridView();
    // * preadapt will mark elements which might vanish due to coarsening
    bool elements_may_be_coarsened = true;
    if (pre_adapt_grid)
      elements_may_be_coarsened = grid_.preAdapt();
    for (auto& data : data_) {
      auto& space = std::get<0>(data).access();
      space.pre_adapt();
    }
    // * each discrete function is associated with persistent storage (keeps local DoF vectors, which can be converted
    //   to DynamicVector<RF>) to keep our data:
    //   - all kept leaf elements might change their indices and
    //   - all coarsened elements might vanish
    // * we also need a container to recall those elements where we need to restrict to
    PersistentContainer<G, bool> restriction_required(grid_, 0, false);
    // * walk the current leaf of the grid ...
    for (auto&& element : elements(grid_view)) {
      for (auto& data : data_) {
        auto& local_function = std::get<3>(data);
        local_function->bind(element);
        //   ... to store the local DoFs ...
        auto& persistent_data = std::get<2>(data);
        persistent_data[element] = XT::LA::convert_to<DynamicVector<RF>>(local_function->dofs());
      }
      //   ... and to mark father elements
      if (element.mightVanish())
        restriction_required[element.father()] = true;
    }
    // * now walk the grid up all coarser levels ...
    if (elements_may_be_coarsened) {
      for (int level = grid_.maxLevel() - 1; level >= 0; --level) {
        auto level_view = grid_.levelGridView(level);
        for (auto&& element : elements(level_view)) {
          // ... to compute restrictions ...
          for (auto& data : data_) {
            const auto& space = std::get<0>(data).access();
            auto& persistent_data = std::get<2>(data);
            if (restriction_required[element])
              space.restrict_to(element, persistent_data);
          }
          // ... and to mark father elements
          if (element.mightVanish())
            restriction_required[element.father()] = true;
        }
      }
    }
  } // ... pre_adapt(...)

  void adapt(const bool adapt_grid = true)
  {
    auto grid_view = grid_.leafGridView();
    if (adapt_grid)
      grid_.adapt();
    // * clean up data structures
    for (auto& data : data_) {
      auto& persistent_data = std::get<2>(data);
      persistent_data.resize();
      persistent_data.shrinkToFit();
    }
    // * update spaces and resize vectors
    for (auto& data : data_) {
      auto& space = std::get<0>(data).access();
      auto& discrete_function = std::get<1>(data).access();
      space.adapt();
      discrete_function.dofs().resize_after_adapt();
    }
    // * get the data back to the discrete function
    for (auto&& element : elements(grid_view)) {
      for (auto& data : data_) {
        auto& space = std::get<0>(data).access();
        const auto& persistent_data = std::get<2>(data);
        auto& local_function = std::get<3>(data);
        local_function->bind(element);
        local_function->dofs().assign_from(space.prolong_onto(element, persistent_data));
      }
    }
  } // ... adapt(...)

  void post_adapt(const bool post_adapt_grid = true, const bool clear = false)
  {
    if (post_adapt_grid)
      grid_.postAdapt();
    for (auto& data : data_) {
      auto& space = std::get<0>(data).access();
      space.post_adapt();
    }
    if (clear)
      data_.clear();
    else {
      auto old_data = std::move(data_);
      data_ = decltype(data_)();
      for (auto& data : old_data)
        this->append(std::get<0>(data).access(), std::get<1>(data).access());
    }
  } // ... post_adapt(...)

private:
  G& grid_;
  std::list<std::tuple<XT::Common::StorageProvider<SpaceType>,
                       XT::Common::StorageProvider<DiscreteFunctionType>,
                       PersistentContainer<G, DynamicVector<RF>>,
                       std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType>>>
      data_;
}; // class AdaptationHelper


template <class V, class GV, size_t r, size_t rC, class RF>
AdaptationHelper<V, GV, r, rC, RF> make_adaptation_helper(typename GV::Grid& grid,
                                                          SpaceInterface<GV, r, rC, RF>& space,
                                                          DiscreteFunction<V, GV, r, rC, RF>& discrete_function)
{
  AdaptationHelper<V, GV, r, rC, RF> helper(grid);
  helper.append(space, discrete_function);
  return helper;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_ADAPTATION_HH
