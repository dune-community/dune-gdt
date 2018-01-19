// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
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


// forwards, needed for the traits
template <class E, class R = double, size_t rT = 1, size_t rCT = 1, size_t rA = rT, size_t rCA = rCT>
class LocalLambdaBinaryVolumeIntegrand;

template <class E, class I, class R = double, size_t rT = 1, size_t rCT = 1, size_t rA = rT, size_t rCA = rCT>
class LocalLambdaBinaryFaceIntegrand;

template <class E, class I, class R = double, size_t rT = 1, size_t rCT = 1, size_t rA = rT, size_t rCA = rCT>
class LocalLambdaQuaternaryFaceIntegrand;

template <class E, class R = double, size_t r = 1, size_t rC = 1>
class LocalLambdaUnaryVolumeIntegrand;


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


template <class E, class I, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
class LocalLambdaBinaryFaceIntegrandTraits
{
  static_assert(XT::Grid::is_entity<E>::value, "");
  static_assert(XT::Grid::is_intersection<I>::value, "");

public:
  typedef LocalLambdaBinaryFaceIntegrand<E, I, R, rT, rCT, rA, rCA> derived_type;
  typedef E EntityType;
  typedef int LocalfunctionTupleType;
  typedef typename EntityType::Geometry::ctype DomainFieldType;
  static const constexpr size_t dimDomain = EntityType::dimension;
};


template <class E, class I, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
class LocalLambdaQuaternaryFaceIntegrandTraits
{
  static_assert(XT::Grid::is_entity<E>::value, "");
  static_assert(XT::Grid::is_intersection<I>::value, "");

public:
  typedef LocalLambdaQuaternaryFaceIntegrand<E, I, R, rT, rCT, rA, rCA> derived_type;
  typedef E EntityType;
  typedef int LocalfunctionTupleType;
  typedef typename EntityType::Geometry::ctype DomainFieldType;
  static const constexpr size_t dimDomain = EntityType::dimension;
};


template <class E, class R, size_t r, size_t rC>
class LocalLambdaUnaryVolumeIntegrandTraits
{
  static_assert(XT::Grid::is_entity<E>::value, "");

public:
  typedef LocalLambdaUnaryVolumeIntegrand<E, R, r, rC> derived_type;
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
}; // class LocalLambdaBinaryVolumeIntegrand


template <class E, class I, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
class LocalLambdaBinaryFaceIntegrand
    : public LocalFaceIntegrandInterface<internal::LocalLambdaBinaryFaceIntegrandTraits<E, I, R, rT, rCT, rA, rCA>, 2>
{
  typedef LocalFaceIntegrandInterface<internal::LocalLambdaBinaryFaceIntegrandTraits<E, I, R, rT, rCT, rA, rCA>, 2>
      BaseType;
  typedef LocalLambdaBinaryFaceIntegrand<E, I, R, rT, rCT, rA, rCA> ThisType;

public:
  typedef internal::LocalLambdaBinaryFaceIntegrandTraits<E, I, R, rT, rCT, rA, rCA> Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::D;
  using BaseType::d;

  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT> TestBaseType;
  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA> AnsatzBaseType;
  typedef FieldVector<D, d - 1> PointType;
  typedef DynamicMatrix<R> LocalMatrixType;

  typedef std::function<size_t(const TestBaseType&, const AnsatzBaseType&)> OrderLambdaType;
  typedef std::function<void(const TestBaseType&, const AnsatzBaseType&, const I&, const PointType&, LocalMatrixType&)>
      EvaluateLambdaType;

  LocalLambdaBinaryFaceIntegrand(OrderLambdaType order_lambda, EvaluateLambdaType evaluate_lambda)
    : order_lambda_(order_lambda)
    , evaluate_lambda_(evaluate_lambda)
  {
  }

  LocalLambdaBinaryFaceIntegrand(const ThisType&) = default;
  LocalLambdaBinaryFaceIntegrand(ThisType&&) = default;

  LocalfunctionTupleType localFunctions(const EntityType& /*entity*/) const
  {
    return 0; // just a dummy
  }

  size_t order(const LocalfunctionTupleType /*local_functions*/,
               const TestBaseType& test_base,
               const AnsatzBaseType& ansatz_base) const
  {
    return order_lambda_(test_base, ansatz_base);
  }

  void evaluate(const LocalfunctionTupleType& /*local_functions*/,
                const TestBaseType& test_base,
                const AnsatzBaseType& ansatz_base,
                const I& intersection,
                const PointType& local_point,
                LocalMatrixType& ret) const
  {
    evaluate_lambda_(test_base, ansatz_base, intersection, local_point, ret);
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
}; // class LocalLambdaBinaryFaceIntegrand


template <class E, class I, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
class LocalLambdaQuaternaryFaceIntegrand
    : public LocalFaceIntegrandInterface<internal::LocalLambdaQuaternaryFaceIntegrandTraits<E, I, R, rT, rCT, rA, rCA>,
                                         4>
{
  typedef LocalFaceIntegrandInterface<internal::LocalLambdaQuaternaryFaceIntegrandTraits<E, I, R, rT, rCT, rA, rCA>, 4>
      BaseType;
  typedef LocalLambdaQuaternaryFaceIntegrand<E, I, R, rT, rCT, rA, rCA> ThisType;

public:
  typedef internal::LocalLambdaQuaternaryFaceIntegrandTraits<E, I, R, rT, rCT, rA, rCA> Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::D;
  using BaseType::d;

  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT> TestBaseType;
  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA> AnsatzBaseType;
  typedef FieldVector<D, d - 1> PointType;
  typedef DynamicMatrix<R> LocalMatrixType;

  typedef std::function<size_t(const TestBaseType&, const AnsatzBaseType&, const TestBaseType&, const AnsatzBaseType&)>
      OrderLambdaType;
  typedef std::function<void(const TestBaseType&,
                             const AnsatzBaseType&,
                             const TestBaseType&,
                             const AnsatzBaseType&,
                             const I&,
                             const PointType&,
                             LocalMatrixType&,
                             LocalMatrixType&,
                             LocalMatrixType&,
                             LocalMatrixType&)>
      EvaluateLambdaType;

  LocalLambdaQuaternaryFaceIntegrand(OrderLambdaType order_lambda, EvaluateLambdaType evaluate_lambda)
    : order_lambda_(order_lambda)
    , evaluate_lambda_(evaluate_lambda)
  {
  }

  LocalLambdaQuaternaryFaceIntegrand(const ThisType&) = default;
  LocalLambdaQuaternaryFaceIntegrand(ThisType&&) = default;

  LocalfunctionTupleType localFunctions(const EntityType& /*entity*/) const
  {
    return 0; // just a dummy
  }

  size_t order(const LocalfunctionTupleType /*local_functions_entity*/,
               const LocalfunctionTupleType /*local_functions_neighbor*/,
               const TestBaseType& test_base_entity,
               const AnsatzBaseType& ansatz_base_entity,
               const TestBaseType& test_base_neighbor,
               const AnsatzBaseType& ansatz_base_neighbor) const
  {
    return order_lambda_(test_base_entity, ansatz_base_entity, test_base_neighbor, ansatz_base_neighbor);
  }

  void evaluate(const LocalfunctionTupleType& /*local_functions_entity*/,
                const LocalfunctionTupleType& /*local_functions_neighbor*/,
                const TestBaseType& test_base_entity,
                const AnsatzBaseType& ansatz_base_entity,
                const TestBaseType& test_base_neighbor,
                const AnsatzBaseType& ansatz_base_neighbor,
                const I& intersection,
                const PointType& local_point,
                LocalMatrixType& ret_entity_entity,
                LocalMatrixType& ret_neighbor_neighbor,
                LocalMatrixType& ret_entity_neighbor,
                LocalMatrixType& ret_neighbor_entity) const
  {
    evaluate_lambda_(test_base_entity,
                     ansatz_base_entity,
                     test_base_neighbor,
                     ansatz_base_neighbor,
                     intersection,
                     local_point,
                     ret_entity_entity,
                     ret_neighbor_neighbor,
                     ret_entity_neighbor,
                     ret_neighbor_entity);
    if (ret_entity_entity.rows() < test_base_entity.size() || ret_entity_entity.cols() < ansatz_base_entity.size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Your evalaute_lambda destroyed ret!\n   "
                     << "ret_entity_entity is expected to be at least of size test_base_entity.size() x "
                        "ansatz_base_entity.size(),\n   "
                     << "do not call ret_entity_entity.resize(...)!\n   "
                     << "test_base_entity.size(): "
                     << test_base_entity.size()
                     << "\n   ansatz_base_entity.size(): "
                     << ansatz_base_entity.size()
                     << "\n   ret_entity_entity.rows(): "
                     << ret_entity_entity.rows()
                     << "\n   ret_entity_entity.cols(): "
                     << ret_entity_entity.cols());
    if (ret_neighbor_neighbor.rows() < test_base_neighbor.size()
        || ret_neighbor_neighbor.cols() < ansatz_base_neighbor.size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Your evalaute_lambda destroyed ret!\n   "
                     << "ret_neighbor_neighbor is expected to be at least of size test_base_neighbor.size() x "
                        "ansatz_base_neighbor.size(),\n   "
                     << "do not call ret_neighbor_neighbor.resize(...)!\n   "
                     << "test_base_neighbor.size(): "
                     << test_base_neighbor.size()
                     << "\n   ansatz_base_neighbor.size(): "
                     << ansatz_base_neighbor.size()
                     << "\n   ret_neighbor_neighbor.rows(): "
                     << ret_neighbor_neighbor.rows()
                     << "\n   ret_neighbor_neighbor.cols(): "
                     << ret_neighbor_neighbor.cols());
    if (ret_entity_neighbor.rows() < test_base_entity.size()
        || ret_entity_neighbor.cols() < ansatz_base_neighbor.size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Your evalaute_lambda destroyed ret!\n   "
                     << "ret_entity_neighbor is expected to be at least of size test_base_entity.size() x "
                        "ansatz_base_neighbor.size(),\n   "
                     << "do not call ret_entity_neighbor.resize(...)!\n   "
                     << "test_base_entity.size(): "
                     << test_base_entity.size()
                     << "\n   ansatz_base_neighbor.size(): "
                     << ansatz_base_neighbor.size()
                     << "\n   ret_entity_neighbor.rows(): "
                     << ret_entity_neighbor.rows()
                     << "\n   ret_entity_neighbor.cols(): "
                     << ret_entity_neighbor.cols());
    if (ret_neighbor_entity.rows() < test_base_neighbor.size()
        || ret_neighbor_entity.cols() < ansatz_base_entity.size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Your evalaute_lambda destroyed ret!\n   "
                     << "ret_neighbor_entity is expected to be at least of size test_base_neighbor.size() x "
                        "ansatz_base_entity.size(),\n   "
                     << "do not call ret_neighbor_entity.resize(...)!\n   "
                     << "test_base_neighbor.size(): "
                     << test_base_neighbor.size()
                     << "\n   ansatz_base_entity.size(): "
                     << ansatz_base_entity.size()
                     << "\n   ret_neighbor_entity.rows(): "
                     << ret_neighbor_entity.rows()
                     << "\n   ret_neighbor_entity.cols(): "
                     << ret_neighbor_entity.cols());
  } // ... evaluate(...)

private:
  const OrderLambdaType order_lambda_;
  const EvaluateLambdaType evaluate_lambda_;
}; // class LocalLambdaQuaternaryFaceIntegrand


template <class E, class R, size_t r, size_t rC>
class LocalLambdaUnaryVolumeIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalLambdaUnaryVolumeIntegrandTraits<E, R, r, rC>, 1>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalLambdaUnaryVolumeIntegrandTraits<E, R, r, rC>, 1> BaseType;
  typedef LocalLambdaUnaryVolumeIntegrand<E, R, r, rC> ThisType;

public:
  typedef internal::LocalLambdaUnaryVolumeIntegrandTraits<E, R, r, rC> Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::D;
  using BaseType::d;

  typedef XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC> TestBaseType;
  typedef FieldVector<D, d> PointType;
  typedef DynamicVector<R> LocalVectorType;

  typedef std::function<size_t(const TestBaseType&)> OrderLambdaType;
  typedef std::function<void(const TestBaseType&, const PointType&, LocalVectorType&)> EvaluateLambdaType;

  LocalLambdaUnaryVolumeIntegrand(OrderLambdaType order_lambda, EvaluateLambdaType evaluate_lambda)
    : order_lambda_(order_lambda)
    , evaluate_lambda_(evaluate_lambda)
  {
  }

  LocalLambdaUnaryVolumeIntegrand(const ThisType&) = default;
  LocalLambdaUnaryVolumeIntegrand(ThisType&&) = default;

  LocalfunctionTupleType localFunctions(const EntityType& /*entity*/) const
  {
    return 0; // just a dummy
  }

  size_t order(const LocalfunctionTupleType& /*local_functions_tuple*/, const TestBaseType& test_base) const
  {
    return order_lambda_(test_base);
  }

  void evaluate(const LocalfunctionTupleType& /*local_functions_tuple*/,
                const TestBaseType& test_base,
                const PointType& local_point,
                LocalVectorType& ret) const
  {
    evaluate_lambda_(test_base, local_point, ret);
    if (ret.size() < test_base.size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Your evalaute_lambda destroyed ret!\n   "
                     << "ret is expected to be at least of size test_base.size(),\n   "
                     << "do not call ret.resize(...)!\n   "
                     << "test_base.size(): "
                     << test_base.size()
                     << "\n   ret.size(): "
                     << ret.size());
  } // ... evaluate(...)

private:
  const OrderLambdaType order_lambda_;
  const EvaluateLambdaType evaluate_lambda_;
}; // class LocalLambdaUnaryVolumeIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LAMBDA_HH
