// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_LOCAL_OPERATORS_LAGRANGE_PROJECTION_HH
#define DUNE_GDT_LOCAL_OPERATORS_LAGRANGE_PROJECTION_HH

#include <type_traits>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/local/discretefunction.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required for the traits
class LocalLagrangeProjectionOperator;


namespace internal {


class LocalLagrangeProjectionOperatorTraits
{
public:
  typedef LocalLagrangeProjectionOperator derived_type;
};


} // namespace internal


class LocalLagrangeProjectionOperator : public LocalOperatorInterface<internal::LocalLagrangeProjectionOperatorTraits>
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
  typedef internal::LocalLagrangeProjectionOperatorTraits Traits;

  /**
   * \brief Applies the Lagrange projection locally.
   */
  template <class E, class D, size_t d, class R, class RangeSpaceType, class VectorType>
  typename std::enable_if<StaticCheck<E, D, d, R, 1, 1, RangeSpaceType, VectorType>::value, void>::type
  apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1, 1>& source,
        LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    const auto& entity = local_range.entity();
    const auto local_source = source.local_function(entity);
    auto& local_DoF_vector = local_range.vector();
    const size_t size = local_DoF_vector.size();
    const auto lagrange_points = local_range.space().lagrange_points(entity);
    assert(lagrange_points.size() == size);
    for (size_t ii = 0; ii < size; ++ii)
      local_DoF_vector.set(ii, local_source->evaluate(lagrange_points[ii]));
  } // ... apply(...)

  /**
   * \brief     Applies the Lagrange projection locally (for vector-valued functions)
   * \attention This is the old implementation which will throw a NotImplemented error. If you need this, test and
   *            correct the implementation!
   */
  template <class E, class D, size_t d, class R, size_t r, class RangeSpaceType, class VectorType>
  typename std::enable_if<StaticCheck<E, D, d, R, r, 1, RangeSpaceType, VectorType>::value, void>::type
  apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, 1>& source,
        LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    DUNE_THROW(NotImplemented, "Think about this!"); // For vector valued functions, the code below has an implicit
    const auto& entity = local_range.entity(); // assumption about the mapper, which may not be ok!
    const auto local_source = source.local_function(entity);
    auto& local_DoF_vector = local_range.vector();
    const auto lagrange_points = local_range.space().lagrange_points(entity);
    // and do the work (see below)
    assert(lagrange_points.size() == local_DoF_vector.size());
    size_t kk = 0;
    for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
      if (std::isinf(local_DoF_vector.get(kk))) { // Assumes that the global DoF vector was set to infinity beforehand,
        // evaluate source function               // which is not the case anymore.
        const auto& lagrange_point = lagrange_points[ii];
        const auto source_value = local_source->evaluate(lagrange_point);
        // set DoFs
        for (size_t jj = 0; jj < r; ++jj, ++kk)
          local_DoF_vector.set(kk, source_value[jj]);
      } else
        kk += r;
    }
  } // ... apply(...)
}; // class LocalLagrangeProjectionOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_LAGRANGE_PROJECTION_HH
