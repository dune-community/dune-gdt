// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_OPERATORS_L1_HH
#define DUNE_GDT_OPERATORS_L1_HH

#include <type_traits>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/gdt/type_traits.hh>

#include "l2.hh"

namespace Dune {
namespace GDT {


// forward, needed for the traits
template <class GridLayer, class Field = double>
class L1Operator;


namespace internal {


template <class GridLayerType, class Field>
class L1OperatorTraits
{
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");

public:
  typedef L1Operator<GridLayerType, Field> derived_type;
  typedef NoJacobian JacobianType;
  typedef Field FieldType;
};


} // namespace internal


template <class GridLayerType, class Field>
class L1Operator : public OperatorInterface<internal::L1OperatorTraits<GridLayerType, Field>>
{
  typedef OperatorInterface<internal::L1OperatorTraits<GridLayerType, Field>> BaseType;

public:
  using typename BaseType::FieldType;
  using typename BaseType::JacobianType;

  using E = XT::Grid::extract_entity_t<GridLayerType>;
  using D = typename GridLayerType::ctype;
  static const size_t d = GridLayerType::dimension;

  L1Operator(GridLayerType grid_layer, const size_t over_integrate = 0)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const Dune::XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class RangeType, class SourceType>
  FieldType
  apply2(const RangeType& range, const SourceType& source, const Dune::XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class SourceType>
  JacobianType jacobian(const SourceType& source, const Dune::XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class SourceType>
  void jacobian(const SourceType& source, JacobianType& jac, const Dune::XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  using BaseType::apply_inverse;

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& range,
                     SourceType& source,
                     const XT::Common::Configuration& opts,
                     const Dune::XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "");
  }

  XT::Common::Configuration invert_options(const std::string& type) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class R, size_t r>
  FieldType induced_norm(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, 1>& range,
                         const Dune::XT::Common::Parameter& param = {}) const
  {
    return std::sqrt(L2Operator<GridLayerType, FieldType>(grid_layer_, over_integrate_)
                         .apply2(XT::Functions::ConstantFunction<E, D, d, R, r>(1.), range, param));
  }

private:
  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class L1Operator


// ///////////////////////// //
// make_weighted_l2_operator //
// ///////////////////////// //

template <class GridLayerType, class R = double>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value, std::unique_ptr<L1Operator<GridLayerType, R>>>::type
make_l1_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L1Operator<GridLayerType, R>>(grid_layer, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_L1_HH
