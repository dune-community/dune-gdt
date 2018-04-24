// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_L2_PROJECTION_HH
#define DUNE_GDT_LOCAL_OPERATORS_L2_PROJECTION_HH

#include <type_traits>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/solver.hh>

#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/spaces/fv/interface.hh>

#include "integrals.hh"
#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required for the traits
class LocalL2ProjectionOperator;


namespace internal {


class LocalL2ProjectionOperatorTraits
{
public:
  typedef LocalL2ProjectionOperator derived_type;
};


} // namespace internal


class LocalL2ProjectionOperator : public LocalOperatorInterface<internal::LocalL2ProjectionOperatorTraits>
{
  template <class E, class D, size_t d, class R, size_t r, size_t rC, class RangeSpaceType, class VectorType>
  struct StaticCheck
  {
    static const bool value = std::is_same<E, typename RangeSpaceType::EntityType>::value
                              && std::is_same<D, typename RangeSpaceType::DomainFieldType>::value
                              && d == RangeSpaceType::dimDomain && r == RangeSpaceType::dimRange
                              && rC == RangeSpaceType::dimRangeCols;
  };

public:
  typedef internal::LocalL2ProjectionOperatorTraits Traits;

  LocalL2ProjectionOperator(const size_t over_integrate = 0, const XT::Common::Parameter& param = {})
    : over_integrate_(over_integrate)
    , param_(param)
  {
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class RangeSpaceType, class VectorType>
  typename std::enable_if<StaticCheck<E, D, d, R, r, rC, RangeSpaceType, VectorType>::value
                              && !is_fv_space<RangeSpaceType>::value,
                          void>::type
  apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
        LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    // create local L2 operator
    typedef XT::Functions::ConstantFunction<E, D, d, R, 1> OneType;
    const OneType one(1.); // <- is not actually used, just needed for the product evaluation
    const LocalVolumeIntegralOperator<LocalProductIntegrand<OneType>,
                                      typename RangeSpaceType::BaseFunctionSetType,
                                      typename RangeSpaceType::BaseFunctionSetType,
                                      R>
        local_l2_operator(over_integrate_, one, param_);
    // and functional
    typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC> SourceType;
    const LocalVolumeIntegralFunctional<LocalProductIntegrand<SourceType>,
                                        typename RangeSpaceType::BaseFunctionSetType,
                                        R>
        local_l2_functional(over_integrate_, source, param_);
    // create local lhs and rhs
    const auto& local_basis = local_range.basis();
    const size_t size = local_basis.size();
    XT::LA::CommonDenseMatrix<R> local_matrix(size, size);
    XT::LA::CommonDenseVector<R> local_vector(size);
    // assemble
    local_l2_operator.apply2(local_basis, local_basis, local_matrix);
    local_l2_functional.apply(local_basis, local_vector.backend());
    // solve
    XT::LA::CommonDenseVector<R> local_solution(size);
    try {
      XT::LA::Solver<XT::LA::CommonDenseMatrix<R>>(local_matrix).apply(local_vector, local_solution);
    } catch (XT::LA::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(projection_error,
                 "L2 projection failed because a local matrix could not be inverted!\n\n"
                     << "This was the original error: "
                     << ee.what());
    }
    // set local DoFs
    auto& local_range_vector = local_range.vector();
    assert(local_range_vector.size() == local_solution.size());
    for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
      local_range_vector.set(ii, local_solution.get_entry(ii));
  } // ... apply(...)

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class RangeSpaceType, class VectorType>
  typename std::enable_if<StaticCheck<E, D, d, R, r, rC, RangeSpaceType, VectorType>::value
                              && is_fv_space<RangeSpaceType>::value,
                          void>::type
  apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
        LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    // create local L2 volume integral functional
    typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC> SourceType;
    const LocalVolumeIntegralFunctional<LocalFVProductIntegrand<SourceType>,
                                        typename RangeSpaceType::BaseFunctionSetType,
                                        R>
        local_l2_functional(over_integrate_, source, param_);
    Dune::DynamicVector<R> local_vector(local_range.basis().size(), 0.);
    const auto& entity = local_range.entity();
    local_l2_functional.apply(local_range.basis(), local_vector);
    local_vector /= entity.geometry().volume();
    auto& local_range_vector = local_range.vector();
    for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
      local_range_vector.set(ii, local_vector[ii]);
  } // ... apply(...) for FV spaces

private:
  const size_t over_integrate_;
  const XT::Common::Parameter param_;
}; // class LocalL2ProjectionOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_L2_PROJECTION_HH
