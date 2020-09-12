// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
#define DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH

#include <memory>
#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/bound-object.hh>

#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


template <class SourceVector,
          class SourceGridView,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceField = double,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeField = SourceField,
          class RangeGridView = SourceGridView,
          class RangeVector = SourceVector>
class LocalElementOperatorInterface
  : public XT::Common::ParametricInterface
  , public XT::Grid::ElementBoundObject<XT::Grid::extract_entity_t<SourceGridView>>
{
  static_assert(
      std::is_same<XT::Grid::extract_entity_t<SourceGridView>, XT::Grid::extract_entity_t<RangeGridView>>::value, "");

  using ThisType = LocalElementOperatorInterface;
  using BaseType = XT::Grid::ElementBoundObject<XT::Grid::extract_entity_t<SourceGridView>>;

public:
  using SV = SourceVector;
  using SGV = SourceGridView;
  static constexpr size_t s_r = source_range_dim;
  static constexpr size_t s_rC = source_range_dim_cols;
  using SR = SourceField;

  using RV = RangeVector;
  using RGV = RangeGridView;
  static constexpr size_t r_r = range_range_dim;
  static constexpr size_t r_rC = range_range_dim_cols;
  using RR = RangeField;
  using LocalRangeType = LocalDiscreteFunction<RV, RGV, r_r, r_rC, RR>;

  static constexpr size_t d = LocalRangeType::d;
  using D = typename LocalRangeType::D;
  using E = typename BaseType::ElementType;

  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, SR>;
  using LocalSourceType = typename SourceType::LocalFunctionType;
  using DiscreteSourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SR>;
  using SourceSpaceType = typename DiscreteSourceType::SpaceType;

  // Allows construction without source, source has to be set by a call to with_source before calling apply
  LocalElementOperatorInterface(const size_t num_local_sources = 1, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , source_()
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < num_local_sources; ++ii)
      local_sources_.emplace_back(nullptr);
  }

  LocalElementOperatorInterface(const SourceType& source,
                                const size_t num_local_sources = 1,
                                const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , source_(source)
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < num_local_sources; ++ii)
      local_sources_.emplace_back(source_.access().local_function());
  }

  LocalElementOperatorInterface(const SourceSpaceType& source_space,
                                const SV& source_vector,
                                const size_t num_local_sources = 1,
                                const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , source_(new DiscreteSourceType(source_space, source_vector))
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < num_local_sources; ++ii)
      local_sources_.emplace_back(source_.access().local_function());
  }

  LocalElementOperatorInterface(const ThisType& other)
    : XT::Common::ParametricInterface(other)
    , XT::Grid::ElementBoundObject<E>(other)
    , source_(other.source_)
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < other.local_sources_.size(); ++ii)
      local_sources_.emplace_back(other.local_sources_[ii] ? source_.access().local_function() : nullptr);
  }

  virtual ~LocalElementOperatorInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  virtual bool linear() const
  {
    return false;
  }

  virtual void apply(LocalRangeType& local_range, const XT::Common::Parameter& param = {}) const = 0;

  virtual std::unique_ptr<ThisType> with_source(const SourceType& src) const
  {
    auto ret = copy();
    ret->source_ = XT::Common::ConstStorageProvider<SourceType>(src);
    for (size_t ii = 0; ii < local_sources_.size(); ++ii)
      ret->local_sources_[ii] = ret->source().local_function();
    return ret;
  }

  const SourceType& source() const
  {
    return source_.access();
  }

  const std::vector<std::unique_ptr<LocalSourceType>>& local_sources() const
  {
    return local_sources_;
  }

protected:
  // We are binding the first local_source to ele and leave the others unbound
  void post_bind(const E& ele) override
  {
    if (local_sources_.size() > 0 && local_sources_[0])
      local_sources_[0]->bind(ele);
  }

  XT::Common::ConstStorageProvider<SourceType> source_;
  std::vector<std::unique_ptr<LocalSourceType>> local_sources_;
}; // class LocalElementOperatorInterface


template <class Intersection,
          class SourceVector,
          class SourceGridView,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceField = double,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeField = SourceField,
          class InsideRangeGridView = SourceGridView,
          class InsideRangeVector = SourceVector,
          class OutsideRangeGridView = InsideRangeGridView,
          class OutsideRangeVector = InsideRangeVector>
class LocalIntersectionOperatorInterface
  : public XT::Common::ParametricInterface
  , public XT::Grid::IntersectionBoundObject<Intersection>
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");
  static_assert(std::is_same<typename Intersection::Entity, XT::Grid::extract_entity_t<InsideRangeGridView>>::value,
                "");
  static_assert(std::is_same<typename Intersection::Entity, XT::Grid::extract_entity_t<OutsideRangeGridView>>::value,
                "");

  using ThisType = LocalIntersectionOperatorInterface;

public:
  static constexpr size_t d = Intersection::Entity::dimension;
  using D = typename Intersection::ctype;
  using I = Intersection;
  using E = typename I::Entity;
  using IntersectionType = Intersection;

  using SV = SourceVector;
  using SGV = SourceGridView;
  static constexpr size_t s_r = source_range_dim;
  static constexpr size_t s_rC = source_range_dim_cols;
  using SF = SourceField;
  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, SF>;
  using LocalSourceType = typename SourceType::LocalFunctionType;
  using DiscreteSourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>;
  using SourceSpaceType = typename DiscreteSourceType::SpaceType;

  using IRV = InsideRangeVector;
  using IRGV = InsideRangeGridView;
  static constexpr size_t r_r = range_range_dim;
  static constexpr size_t r_rC = range_range_dim_cols;
  using RF = RangeField;
  using LocalInsideRangeType = LocalDiscreteFunction<IRV, IRGV, r_r, r_rC, RF>;

  using ORV = OutsideRangeVector;
  using ORGV = OutsideRangeGridView;
  using LocalOutsideRangeType = LocalDiscreteFunction<ORV, ORGV, r_r, r_rC, RF>;

  // Allows construction without source, source has to be set by a call to with_source before calling apply
  LocalIntersectionOperatorInterface(const size_t num_local_sources = 2,
                                     const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , source_()
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < num_local_sources; ++ii)
      local_sources_.emplace_back(nullptr);
  }

  LocalIntersectionOperatorInterface(const SourceType& src,
                                     const size_t num_local_sources = 2,
                                     const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , source_(src)
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < num_local_sources; ++ii)
      local_sources_.emplace_back(source_.access().local_function());
  }

  LocalIntersectionOperatorInterface(const SourceSpaceType& source_space,
                                     const SV& source_vector,
                                     const size_t num_local_sources = 2,
                                     const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , source_(new DiscreteSourceType(source_space, source_vector))
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < num_local_sources; ++ii)
      local_sources_.emplace_back(source_.access().local_function());
  }

  LocalIntersectionOperatorInterface(const ThisType& other)
    : XT::Common::ParametricInterface(other)
    , XT::Grid::IntersectionBoundObject<Intersection>(other)
    , source_(other.source_)
    , local_sources_(0)
  {
    for (size_t ii = 0; ii < other.local_sources_.size(); ++ii)
      local_sources_.emplace_back(other.local_sources_[ii] ? source_.access().local_function() : nullptr);
  }

  virtual ~LocalIntersectionOperatorInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  virtual bool linear() const
  {
    return false;
  }

  /**
   * \note Presumes that local_range_inside is already bound to intersection.inside() and local_range_outside is
   *       already bound to intersection.outside()!
   **/
  virtual void apply(LocalInsideRangeType& local_range_inside,
                     LocalOutsideRangeType& local_range_outside,
                     const XT::Common::Parameter& param = {}) const = 0;

  virtual std::unique_ptr<ThisType> with_source(const SourceType& src) const
  {
    auto ret = copy();
    ret->source_ = XT::Common::ConstStorageProvider<SourceType>(src);
    for (size_t ii = 0; ii < local_sources_.size(); ++ii)
      ret->local_sources_[ii] = ret->source().local_function();
    return ret;
  }

  const SourceType& source() const
  {
    return source_.access();
  }

  const std::vector<std::unique_ptr<LocalSourceType>>& local_sources() const
  {
    return local_sources_;
  }

protected:
  // We are binding the first local_source to intersection.inside() and the second one to intersection.outside()
  void post_bind(const I& inter) override
  {
    if (local_sources_.size() > 0 && local_sources_[0])
      local_sources_[0]->bind(inter.inside());
    if (local_sources_.size() > 1 && local_sources_[1])
      local_sources_[1]->bind(inter.neighbor() ? inter.outside() : inter.inside());
  }

  XT::Common::ConstStorageProvider<SourceType> source_;
  std::vector<std::unique_ptr<LocalSourceType>> local_sources_;
}; // class LocalIntersectionOperatorInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
