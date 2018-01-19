// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_FV_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_FV_HH

#include <tuple>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsIntegrand;

template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsJacobianIntegrand;

template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsNewtonIntegrand;

namespace internal {


/**
 * \tparam RhsEvaluationImp Right-hand side function derived from RhsEvaluationInterface
 * \tparam SourceType DiscreteFunction that contains the values u_i
 */
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsIntegrandTraits
{
public:
  typedef RhsEvaluationImp RhsEvaluationType;
  typedef typename RhsEvaluationImp::EntityType EntityType;
  typedef typename RhsEvaluationImp::DomainFieldType DomainFieldType;
  static const size_t dimDomain = RhsEvaluationImp::dimDomain;
  typedef LocalFvRhsIntegrand<RhsEvaluationImp, SourceType> derived_type;
  typedef std::tuple<std::unique_ptr<typename RhsEvaluationType::LocalfunctionType>,
                     std::unique_ptr<typename SourceType::LocalfunctionType>,
                     typename SourceType::SpaceType::GridLayerType::ctype>
      LocalfunctionTupleType;
}; // class LocalFvRhsIntegrandTraits


/**
 * \tparam RhsEvaluationImp Right-hand side function derived from RhsEvaluationInterface
 * \tparam SourceType DiscreteFunction that contains the values u_i
 */
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsJacobianIntegrandTraits
{
public:
  typedef RhsEvaluationImp RhsEvaluationType;
  typedef typename RhsEvaluationImp::EntityType EntityType;
  typedef typename RhsEvaluationImp::DomainFieldType DomainFieldType;
  static const size_t dimDomain = RhsEvaluationImp::dimDomain;
  typedef LocalFvRhsJacobianIntegrand<RhsEvaluationImp, SourceType> derived_type;
  typedef std::tuple<std::unique_ptr<typename RhsEvaluationType::LocalfunctionType>,
                     std::unique_ptr<typename SourceType::LocalfunctionType>,
                     typename SourceType::SpaceType::GridLayerType::ctype>
      LocalfunctionTupleType;
}; // class LocalFvRhsJacobianIntegrandTraits

/**
 * \tparam RhsEvaluationImp Right-hand side function derived from RhsEvaluationInterface
 * \tparam SourceType DiscreteFunction that contains the values u_i
 */
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsNewtonIntegrandTraits
{
public:
  typedef RhsEvaluationImp RhsEvaluationType;
  typedef typename RhsEvaluationImp::EntityType EntityType;
  typedef typename RhsEvaluationImp::DomainFieldType DomainFieldType;
  static const size_t dimDomain = RhsEvaluationImp::dimDomain;
  typedef LocalFvRhsNewtonIntegrand<RhsEvaluationImp, SourceType> derived_type;
  typedef std::tuple<size_t,
                     std::unique_ptr<typename SourceType::LocalfunctionType>,
                     typename SourceType::SpaceType::GridLayerType::ctype>
      LocalfunctionTupleType;
}; // class LocalFvRhsNewtonIntegrandTraits


} // namespace internal


/**
 * \tparam RhsEvaluationImp Right-hand side function derived from RhsEvaluationInterface
 * \tparam SourceType DiscreteFunction that contains the values u_i
 */
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalFvRhsIntegrandTraits<RhsEvaluationImp, SourceType>, 1>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalFvRhsIntegrandTraits<RhsEvaluationImp, SourceType>, 1> BaseType;
  typedef LocalFvRhsIntegrand<RhsEvaluationImp, SourceType> ThisType;

public:
  typedef internal::LocalFvRhsIntegrandTraits<RhsEvaluationImp, SourceType> Traits;
  typedef typename Traits::RhsEvaluationType RhsEvaluationType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename RhsEvaluationType::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;

public:
  LocalFvRhsIntegrand(const RhsEvaluationType& rhs_evaluation, const SourceType& source)
    : rhs_evaluation_(rhs_evaluation)
    , source_(source)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(
        rhs_evaluation_.local_function(entity), source_.local_function(entity), entity.geometry().volume());
  }

  template <class R, size_t r, size_t rC>
  size_t
  order(const LocalfunctionTupleType& local_functions_tuple,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& /*test_base*/)
      const
  {
    return std::get<0>(local_functions_tuple)->order({});
  }

  template <class R, size_t r, size_t rC>
  void evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& /*test_base*/,
      const Dune::FieldVector<DomainFieldType, dimDomain>& x_local,
      Dune::DynamicVector<R>& ret) const
  {
    const auto u = std::get<1>(local_functions_tuple)->evaluate(x_local);
    ret = DynamicVector<R>(std::get<0>(local_functions_tuple)->evaluate(x_local, u));
    ret /= std::get<2>(local_functions_tuple);
  }

private:
  const RhsEvaluationType& rhs_evaluation_;
  const SourceType& source_;
}; // class LocalFvRhsIntegrand


/**
 * \tparam RhsEvaluationImp Right-hand side function derived from RhsEvaluationInterface
 * \tparam SourceType DiscreteFunction that contains the values u_i
 */
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsJacobianIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalFvRhsJacobianIntegrandTraits<RhsEvaluationImp, SourceType>, 2>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalFvRhsJacobianIntegrandTraits<RhsEvaluationImp, SourceType>, 2>
      BaseType;
  typedef LocalFvRhsJacobianIntegrand<RhsEvaluationImp, SourceType> ThisType;

public:
  typedef internal::LocalFvRhsJacobianIntegrandTraits<RhsEvaluationImp, SourceType> Traits;
  typedef typename Traits::RhsEvaluationType RhsEvaluationType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename RhsEvaluationType::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;

public:
  LocalFvRhsJacobianIntegrand(const RhsEvaluationType& rhs_evaluation, const SourceType& source)
    : rhs_evaluation_(rhs_evaluation)
    , source_(source)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(
        rhs_evaluation_.local_function(entity), source_.local_function(entity), entity.geometry().volume());
  }

  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
      const XT::Functions::
          LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/) const
  {
    return std::get<0>(local_functions_tuple)->order({});
  }

  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*test_base*/,
      const XT::Functions::
          LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/,
      const Dune::FieldVector<DomainFieldType, dimDomain>& x_local,
      Dune::DynamicMatrix<R>& ret) const
  {
    const auto u = std::get<1>(local_functions_tuple)->evaluate(x_local);
    ret = std::get<0>(local_functions_tuple)->partial_u(x_local, u);
    ret /= std::get<2>(local_functions_tuple);
  }

private:
  const RhsEvaluationType& rhs_evaluation_;
  const SourceType& source_;
}; // class LocalFvRhsJacobianIntegrand

/**
 * \tparam RhsEvaluationImp Right-hand side function derived from RhsEvaluationInterface
 * \tparam SourceType DiscreteFunction that contains the values u_i
 */
template <class RhsEvaluationImp, class SourceType>
class LocalFvRhsNewtonIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalFvRhsNewtonIntegrandTraits<RhsEvaluationImp, SourceType>, 2>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalFvRhsNewtonIntegrandTraits<RhsEvaluationImp, SourceType>, 2>
      BaseType;
  typedef LocalFvRhsNewtonIntegrand<RhsEvaluationImp, SourceType> ThisType;

public:
  typedef internal::LocalFvRhsNewtonIntegrandTraits<RhsEvaluationImp, SourceType> Traits;
  typedef typename Traits::RhsEvaluationType RhsEvaluationType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename RhsEvaluationType::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;

public:
  LocalFvRhsNewtonIntegrand(const RhsEvaluationType& rhs_evaluation,
                            const SourceType& source,
                            const XT::Common::Parameter param)
    : rhs_evaluation_(rhs_evaluation)
    , source_(source)
    , factor_(param.get("factor")[0])
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(rhs_evaluation_.order(entity), source_.local_function(entity), entity.geometry().volume());
  }

  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
      const XT::Functions::
          LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/) const
  {
    return std::get<0>(local_functions_tuple);
  }

  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
      const XT::Functions::
          LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/,
      const Dune::FieldVector<DomainFieldType, dimDomain>& x_local,
      Dune::DynamicMatrix<R>& ret) const
  {
    const auto& entity = test_base.entity();
    const auto u = std::get<1>(local_functions_tuple)->evaluate(x_local);
    ret = rhs_evaluation_.jacobian(u, entity, x_local);
    ret *= factor_ / std::get<2>(local_functions_tuple);
    for (size_t ii = 0; ii < ret.size(); ++ii)
      ret[ii][ii] += 1.;
  }

private:
  const RhsEvaluationType& rhs_evaluation_;
  const SourceType& source_;
  const double factor_;
}; // class LocalFvRhsNewtonIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_FV_HH
