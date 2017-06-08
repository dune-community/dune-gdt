// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_LAMBDA_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_LAMBDA_HH

#include <functional>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, needed for the traits
template <class E, class R = double, size_t rT = 1, size_t rCT = 1, size_t rA = rT, size_t rCA = rCT>
class LocalLambdaBinaryVolumeIntegrand;


namespace internal {


template <class E, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
class LocalLambdaBinaryVolumeIntegrandTraits
{
  static_assert(XT::Grid::is_entity<E>::value, "");

public:
  typedef LocalLambdaBinaryVolumeIntegrand<E, R, rT, rCT, rA, rCA> derived_type;
  typedef E EntityType;
  typedef int LocalfunctionTupleType;
  typedef typename EntityType::Geometry::ctype DomainFieldType;
  static const constexpr size_t dimDomain = EntityType::dimension;
};


} // namespace internal


template <class E, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
class LocalLambdaBinaryVolumeIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalLambdaBinaryVolumeIntegrandTraits<E, R, rT, rCT, rA, rCA>, 2>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalLambdaBinaryVolumeIntegrandTraits<E, R, rT, rCT, rA, rCA>, 2>
      BaseType;
  typedef LocalLambdaBinaryVolumeIntegrand<E, R, rT, rCT, rA, rCA> ThisType;

public:
  typedef internal::LocalLambdaBinaryVolumeIntegrandTraits<E, R, rT, rCT, rA, rCA> Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::D;
  using BaseType::d;

  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT> TestBaseType;
  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA> AnsatzBaseType;
  typedef FieldVector<D, d> PointType;
  typedef DynamicMatrix<R> LocalMatrixType;

  typedef std::function<size_t(const TestBaseType&, const AnsatzBaseType&)> OrderLambdaType;
  typedef std::function<void(const TestBaseType&, const AnsatzBaseType&, const PointType&, LocalMatrixType&)>
      EvaluateLambdaType;

  LocalLambdaBinaryVolumeIntegrand(OrderLambdaType order_lambda, EvaluateLambdaType evaluate_lambda)
    : order_lambda_(order_lambda)
    , evaluate_lambda_(evaluate_lambda)
  {
  }

  LocalLambdaBinaryVolumeIntegrand(const ThisType&) = default;
  LocalLambdaBinaryVolumeIntegrand(ThisType&&) = default;

  LocalfunctionTupleType localFunctions(const EntityType& /*entity*/) const
  {
    return 0; // just a dummy
  }

  size_t order(const LocalfunctionTupleType& /*local_functions_tuple*/,
               const TestBaseType& test_base,
               const AnsatzBaseType& ansatz_base) const
  {
    return order_lambda_(test_base, ansatz_base);
  }

  void evaluate(const LocalfunctionTupleType& /*local_functions_tuple*/,
                const TestBaseType& test_base,
                const AnsatzBaseType& ansatz_base,
                const PointType& local_point,
                LocalMatrixType& ret) const
  {
    evaluate_lambda_(test_base, ansatz_base, local_point, ret);
    if (ret.rows() < test_base.size() || ret.cols() < ansatz_base.size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Your evalaute_lambda destroyed ret!\n   "
                     << "ret is expected to be at least of size test_base.size() x ansatz_base.size(),\n   "
                     << "do not call ret.resize(...)!\n   "
                     << "test_base.size(): "
                     << test_base.size()
                     << "\n   ansatz_base.size(): "
                     << ansatz_base.size()
                     << "\n   ret.rows(): "
                     << ret.rows()
                     << "\n   ret.cols(): "
                     << ret.cols());
  } // ... evaluate(...)

private:
  const OrderLambdaType order_lambda_;
  const EvaluateLambdaType evaluate_lambda_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LAMBDA_HH
