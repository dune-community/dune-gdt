// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_ADAPTATION_HH
#define DUNE_GDT_DISCRETEFUNCTION_ADAPTATION_HH

#include <list>

#include <dune/common/dynvector.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/timedlogging.hh>
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
class AdaptationHelper : public XT::Common::WithLogger<AdaptationHelper<V, GV, r, rC, RF>>
{
  using ThisType = AdaptationHelper;
  using Logger = XT::Common::WithLogger<AdaptationHelper<V, GV, r, rC, RF>>;

public:
  using DiscreteFunctionType = DiscreteFunction<V, GV, r, rC, RF>;
  using SpaceType = SpaceInterface<GV, r, rC, RF>;
  using G = typename GV::Grid;
  static_assert(!XT::Grid::is_yaspgrid<G>::value, "The PersistentContainer is known to segfault for YaspGrid!");

  AdaptationHelper(G& grd, const std::string& logging_prefix = "")
    : Logger(logging_prefix.empty() ? "gdt" : "gdt.tools.adaptation-helper",
             logging_prefix.empty() ? "AdaptationHelper" : logging_prefix,
             /*logging_disabled=*/logging_prefix.empty())
    , grid_(grd)
    , data_(new std::remove_reference_t<decltype(*data_)>)
  {
    LOG_(info) << this->logging_id << "(&grd=" << &grd << ")" << std::endl;
  }

  ThisType& append(SpaceType& space, DiscreteFunctionType& discrete_function)
  {
    LOG_(info) << this->logging_id << ".append(space=" << space << ", discrete_function=" << &discrete_function << ")"
               << std::endl;
    data_->emplace_back(XT::Common::StorageProvider<SpaceType>(space),
                        XT::Common::StorageProvider<DiscreteFunctionType>(discrete_function),
                        PersistentContainer<G, std::pair<DynamicVector<RF>, DynamicVector<RF>>>(
                            grid_, 0, std::make_pair(DynamicVector<RF>(), DynamicVector<RF>())),
                        std::move(discrete_function.local_discrete_function()));
    return *this;
  }

  void pre_adapt(const bool pre_adapt_grid = true)
  {
    LOG_(info) << this->logging_id << ".pre_adapt(pre_adapt_grid=" << pre_adapt_grid << ")" << std::endl;
    auto grid_view = grid_.leafGridView();
    // * preadapt will mark elements which might vanish due to coarsening
    bool elements_may_be_coarsened = true;
    if (pre_adapt_grid) {
      LOG_(info) << "    pre-adapting grid ..." << std::endl;
      elements_may_be_coarsened = grid_.preAdapt();
    }
    LOG_(info) << "    pre-adapting " << data_->size() << " spaces ..." << std::endl;
    for (auto& data : *data_) {
      auto& space = std::get<0>(data).access();
      space.pre_adapt();
    }
    LOG_(info) << "    storing persistent leaf data ..." << std::endl;
    // * each discrete function is associated with persistent storage (see data_, keeps local DoF vectors, which can be
    //   converted to DynamicVector<RF>) to keep our data:
    //   - all kept leaf elements might change their indices and
    //   - all coarsened elements might vanish
    // * we also need a container to recall those elements where we need to restrict to
    //   [note: we would like to use `PersistentContainer<G, bool> restriction_required(grid_, 0, false);` here, but
    //    then `restriction_required[element.father()] = true` does not compile for some grids]
    PersistentContainer<G, int> restriction_required(grid_, 0, 0);
    // * walk the current leaf of the grid ...
    for (auto&& element : elements(grid_view)) {
      for (auto& data : *data_) {
        //   ... to get the local FE data ...
        auto local_FE_data = std::get<0>(data).access().basis().localize(element)->backup();
        //   ... and the local DoF data ...
        auto& local_function = std::get<3>(data);
        local_function->bind(element);
        auto local_DoF_data = XT::LA::convert_to<DynamicVector<RF>>(local_function->dofs());
        //   ... and to store them ...
        auto& persistent_data = std::get<2>(data);
        persistent_data[element] = std::make_pair(std::move(local_FE_data), std::move(local_DoF_data));
      }
      //   ... and to mark father elements
      if (element.mightVanish())
        restriction_required[element.father()] = 1;
    }
    LOG_(info) << "    computing restrictions ..." << std::endl;
    // * now walk the grid up all coarser levels ...
    if (elements_may_be_coarsened) {
      for (int level = grid_.maxLevel() - 1; level >= 0; --level) {
        auto level_view = grid_.levelGridView(level);
        for (auto&& element : elements(level_view)) {
          // ... to compute restrictions ...
          for (auto& data : *data_) {
            const auto& space = std::get<0>(data).access();
            auto& persistent_data = std::get<2>(data);
            if (restriction_required[element])
              space.restrict_to(element, persistent_data);
          }
          // ... and to mark father elements
          if (element.mightVanish()) {
            DUNE_THROW_IF(
                level == 0, Exceptions::space_error, "It does not make sense that a level 0 element might vanish!!");
            restriction_required[element.father()] = true;
          }
        }
      }
    }
  } // ... pre_adapt(...)

  void adapt(const bool adapt_grid = true)
  {
    LOG_(info) << this->logging_id << ".adapt(adapt_grid=" << adapt_grid << ")" << std::endl;
    auto grid_view = grid_.leafGridView();
    if (adapt_grid) {
      LOG_(info) << "    adapting grid ..." << std::endl;
      grid_.adapt();
    }
    LOG_(info) << "    adapting persistent data ..." << std::endl;
    // * clean up data structures
    for (auto& data : *data_) {
      auto& persistent_data = std::get<2>(data);
      persistent_data.resize();
      persistent_data.shrinkToFit();
    }
    LOG_(info) << "    adapting " << data_->size() << " spaces ..." << std::endl;
    // * update spaces and resize vectors
    for (auto& data : *data_) {
      auto& space = std::get<0>(data).access();
      auto& discrete_function = std::get<1>(data).access();
      space.adapt();
      discrete_function.dofs().resize_after_adapt();
    }
    LOG_(info) << "    computing prolongations ..." << std::endl;
    // * get the data back to the discrete function
    for (auto&& element : elements(grid_view)) {
      for (auto& data : *data_) {
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
    LOG_(info) << this->logging_id << ".post_adapt(post_adapt_grid=" << post_adapt_grid << ", clear=" << clear << ")"
               << std::endl;
    if (post_adapt_grid) {
      LOG_(info) << "    post-adapting grid ..." << std::endl;
      grid_.postAdapt();
    }
    LOG_(info) << "    post-adapting " << data_->size() << " spaces ..." << std::endl;
    for (auto& data : *data_) {
      auto& space = std::get<0>(data).access();
      space.post_adapt();
    }
    if (clear) {
      LOG_(info) << "    clearing data ..." << std::endl;
      data_->clear();
    } else {
      LOG_(info) << "    keeping track of " << data_->size() << " spaces:" << std::endl;
      auto old_data = data_;
      data_ = std::make_shared<std::remove_reference_t<decltype(*data_)>>();
      for (auto& data : *old_data) {
        this->append(std::get<0>(data).access(), std::get<1>(data).access());
      }
    }
  } // ... post_adapt(...)

protected:
  G& grid_;
  std::shared_ptr<std::list<std::tuple<XT::Common::StorageProvider<SpaceType>,
                                       XT::Common::StorageProvider<DiscreteFunctionType>,
                                       PersistentContainer<G, std::pair<DynamicVector<RF>, DynamicVector<RF>>>,
                                       std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType>>>>
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
