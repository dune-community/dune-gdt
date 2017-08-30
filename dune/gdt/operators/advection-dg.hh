// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_OPERATORS_ADVECTION_DG_HH

#include <functional>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/local/operators/advection-dg.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// forward
template <class DF, class GL>
class AdvectionDgOperator;


namespace internal {


template <class DF, class GL>
class AdvectionDgOperatorTraits
{
public:
  using derived_type = AdvectionDgOperator<DF, GL>;
  using FieldType = double;
  using JacobianType = int;
};


} // namespace internal


template <class DF, class GL = typename DF::SpaceType::GridLayerType>
class AdvectionDgOperator : public OperatorInterface<internal::AdvectionDgOperatorTraits<DF, GL>>
{
  static_assert(is_discrete_function<DF>::value, "");
  using SpaceType = typename DF::SpaceType;
  static_assert(XT::Grid::is_layer<GL>::value, "");
  using ThisType = AdvectionDgOperator<DF, GL>;
  using BaseType = OperatorInterface<internal::AdvectionDgOperatorTraits<DF, GL>>;
  using LocalOperatorType = LocalAdvectionDgOperator<SpaceType, GL>;

public:
  using typename BaseType::FieldType;
  using typename BaseType::JacobianType;
  using FluxType = typename LocalOperatorType::FluxType;
  using NumericalCouplingFluxType = typename LocalOperatorType::NumericalFluxType;

  AdvectionDgOperator(const GL& grid_layer,
                      const FluxType& flux,
                      NumericalCouplingFluxType numerical_coupling_flux = /*simple upwinding*/
                      [](const auto& f, const auto& u, const auto& v, const auto& n) {
                        const auto df = f.partial_u({}, (u + v) / 2.);
                        if ((df[0] * n) > 0)
                          return f.evaluate({}, u) * n;
                        else
                          return f.evaluate({}, v) * n;
                      })
    : grid_layer_(grid_layer)
    , local_operator_(grid_layer_, flux, numerical_coupling_flux)
  {
  }

  AdvectionDgOperator(const ThisType& other) = delete;
  AdvectionDgOperator(ThisType&& source) = delete;

  void apply(const DF& source, DF& range, const XT::Common::Parameter& /*param*/ = {}) const
  {
    range.vector() *= 0.;
    LocalizableOperatorBase<GL, DF, DF> walker(grid_layer_, source, range);
    walker.append(local_operator_);
    walker.walk();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/,
                   const SourceType& /*source*/,
                   const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
    return FieldType();
  }

  template <class SourceType>
  JacobianType jacobian(const SourceType& /*source*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
    return JacobianType();
  }

  template <class SourceType>
  void
  jacobian(const SourceType& /*source*/, JacobianType& /*jac*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/,
                     SourceType& /*source*/,
                     const XT::Common::Configuration& /*opts*/,
                     const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "");
    return std::vector<std::string>();
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "");
    return XT::Common::Configuration();
  }

private:
  const GL& grid_layer_;
  const LocalOperatorType local_operator_;
}; // class AdvectionDgOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_DG_HH
