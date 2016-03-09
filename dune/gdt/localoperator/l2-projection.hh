// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_L2_PROJECTION_HH
#define DUNE_GDT_LOCALOPERATOR_L2_PROJECTION_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/discretefunction/local.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localfunctional/integrals.hh>
#include <dune/gdt/localoperator/integrals.hh>
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


class LocalL2ProjectionOperator
  : public LocalOperatorInterface< internal::LocalL2ProjectionOperatorTraits >
{
  template< class E, class D, size_t d, class R, size_t r, size_t rC, class RangeSpaceType, class VectorType >
  struct StaticCheck
  {
    static const bool value = std::is_same< E, typename RangeSpaceType::EntityType >::value
                           && std::is_same< D, typename RangeSpaceType::DomainFieldType >::value
                           && d  == RangeSpaceType::dimDomain
                           && r  == RangeSpaceType::dimRange
                           && rC == RangeSpaceType::dimRangeCols;
  };

public:
  typedef internal::LocalL2ProjectionOperatorTraits Traits;

  LocalL2ProjectionOperator(const size_t over_integrate = 0)
    : over_integrate_(over_integrate)
  {}

  template< class E, class D, size_t d, class R, size_t r, size_t rC, class RangeSpaceType, class VectorType >
      typename std::enable_if< StaticCheck< E, D, d, R, r, rC, RangeSpaceType, VectorType >::value && !(is_fv_space< RangeSpaceType >::value || is_product_fv_space< RangeSpaceType >::value), void >::type
  apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
        LocalDiscreteFunction< RangeSpaceType, VectorType >& local_range) const
  {
    // create local L2 operator
    typedef Stuff::Functions::Constant< E, D, d, R, 1 > OneType;
    const OneType one(1.); // <- is not actually used, just needed for the product evaluation
    const  LocalVolumeIntegralOperator< LocalEvaluation::Product< OneType > >
        local_l2_operator(over_integrate_, one);
    // and functional
    typedef Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC > SourceType;
    const LocalVolumeIntegralFunctional< LocalEvaluation::Product< SourceType > >
        local_l2_functional(over_integrate_, source);
    // create local lhs and rhs
    const auto& local_basis = local_range.basis();
    const size_t size = local_basis.size();
    Stuff::LA::CommonDenseMatrix< R > local_matrix(size, size);
    Stuff::LA::CommonDenseVector< R > local_vector(size);
    // assemble
    local_l2_operator.apply2(local_basis,  local_basis, local_matrix.backend());
    local_l2_functional.apply(local_basis, local_vector.backend());
    // solve
    Stuff::LA::CommonDenseVector< R > local_solution(size);
    try {
      Stuff::LA::Solver< Stuff::LA::CommonDenseMatrix< R > >(local_matrix).apply(local_vector, local_solution);
    } catch (Stuff::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(Exceptions::projection_error,
                 "L2 projection failed because a local matrix could not be inverted!\n\n"
                 << "This was the original error: " << ee.what());
    }
    // set local DoFs
    auto& local_range_vector = local_range.vector();
    assert(local_range_vector.size() == local_solution.size());
    for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
      local_range_vector.set(ii, local_solution[ii]);
  } // ... apply(...)

  // TODO: do not use product evaluation to avoid a lot of multiplications with 0
  template< class E, class D, size_t d, class R, size_t r, size_t rC, class RangeSpaceType, class VectorType >
  typename std::enable_if< StaticCheck< E, D, d, R, r, rC, RangeSpaceType, VectorType >::value && (is_fv_space< RangeSpaceType >::value || is_product_fv_space< RangeSpaceType >::value), void >::type
  apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
        LocalDiscreteFunction< RangeSpaceType, VectorType >& local_range) const
  {
    // create local L2 volume integral functional
    typedef Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC > SourceType;
    const LocalVolumeIntegralFunctional< LocalEvaluation::Product< SourceType > >
        local_l2_functional(over_integrate_, source);
    Stuff::LA::CommonDenseVector< R > local_vector(local_range.basis().size());
    const auto& entity = local_range.entity();
    local_l2_functional.apply(local_range.basis(), local_vector.backend());
    local_vector /= entity.geometry().volume();
    auto& local_range_vector = local_range.vector();
    for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
      local_range_vector.set(ii, local_vector[ii]);
  } // ... apply(...) for FV spaces

private:
  const size_t over_integrate_;
}; // class LocalL2ProjectionOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_L2_PROJECTION_HH
