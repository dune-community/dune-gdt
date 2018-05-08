// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_FV_HH
#define DUNE_GDT_LOCAL_OPERATORS_FV_HH

#include <dune/common/fvector.hh>

#include <dune/xt/common/fvector.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class NumericalFluxType>
class LocalCouplingFvOperator;

template <class NumericalFluxType>
class LocalBoundaryFvOperator;


namespace internal {


template <class NumericalFluxType>
struct LocalCouplingFvOperatorTraits
{
  typedef LocalCouplingFvOperator<NumericalFluxType> derived_type;
};

template <class NumericalFluxType>
struct LocalBoundaryFvOperatorTraits
{
  typedef LocalBoundaryFvOperator<NumericalFluxType> derived_type;
};


} // namespace internal


template <class NumericalFluxType>
class LocalCouplingFvOperator
    : public LocalCouplingOperatorInterface<internal::LocalCouplingFvOperatorTraits<NumericalFluxType>>
{
public:
  typedef typename NumericalFluxType::RangeFieldType RangeFieldType;
  static const size_t dimDomain = NumericalFluxType::dimDomain;
  static const size_t dimRange = NumericalFluxType::dimRange;
  typedef Dune::QuadratureRule<RangeFieldType, dimDomain - 1> QuadratureType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> ResultType;


  template <class... Args>
  explicit LocalCouplingFvOperator(const QuadratureType& quadrature, Args&&... args)
    : quadrature_(quadrature)
    , numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor,
             const XT::Common::Parameter& /*mu*/ = {}) const
  {
    const auto entity = intersection.inside();
    const auto neighbor = intersection.outside();
    const auto local_source_entity = source.local_function(entity);
    const auto local_source_neighbor = source.local_function(neighbor);
    const auto local_functions_tuple_entity = numerical_flux_.local_functions(entity);
    const auto local_functions_tuple_neighbor = numerical_flux_.local_functions(neighbor);
    ResultType result(0);
    for (const auto& quad_point : quadrature_) {
      result += ResultType(numerical_flux_.evaluate(local_functions_tuple_entity,
                                                    local_functions_tuple_neighbor,
                                                    *local_source_entity,
                                                    *local_source_neighbor,
                                                    intersection,
                                                    quad_point.position()))
                * quad_point.weight() * intersection.geometry().integrationElement(quad_point.position());
    }
    local_range_entity.vector().add(result * (1.0 / entity.geometry().volume()));
    local_range_neighbor.vector().add(result * (-1.0 / neighbor.geometry().volume()));
  }

private:
  const QuadratureType& quadrature_;
  const NumericalFluxType numerical_flux_;
};


template <class NumericalFluxType>
class LocalBoundaryFvOperator
    : public LocalBoundaryOperatorInterface<internal::LocalBoundaryFvOperatorTraits<NumericalFluxType>>
{
public:
  typedef typename NumericalFluxType::RangeFieldType RangeFieldType;
  static const size_t dimDomain = NumericalFluxType::dimDomain;
  static const size_t dimRange = NumericalFluxType::dimRange;
  typedef Dune::QuadratureRule<RangeFieldType, dimDomain - 1> QuadratureType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> ResultType;

  template <class... Args>
  explicit LocalBoundaryFvOperator(const QuadratureType& quadrature, Args&&... args)
    : quadrature_(quadrature)
    , numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity) const
  {
    const auto entity = intersection.inside();
    const auto local_source_entity = source.local_function(entity);
    const auto local_functions_tuple = numerical_flux_.local_functions(entity);
    ResultType result(0);
    for (const auto& quad_point : quadrature_) {
      result += ResultType(numerical_flux_.evaluate(
                    local_functions_tuple, *local_source_entity, intersection, quad_point.position()))
                * quad_point.weight() * intersection.geometry().integrationElement(quad_point.position());
    }
    result /= entity.geometry().volume();
    local_range_entity.vector().add(result);
  }

private:
  const QuadratureType& quadrature_;
  const NumericalFluxType numerical_flux_;
}; // class LocalBoundaryFvOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_FV_HH
