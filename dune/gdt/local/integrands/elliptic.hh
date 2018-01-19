// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH

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
template <class DiffusionFactorImp, class DiffusionTensorImp = void>
class LocalEllipticIntegrand;


namespace internal {


/**
 * \brief Traits for the Elliptic evaluation (variant for given diffusion factor and tensor).
 */
template <class DiffusionFactorImp, class DiffusionTensorImp>
class LocalEllipticIntegrandTraits
{
  static_assert(XT::Functions::is_localizable_function<DiffusionFactorImp>::value,
                "DiffusionFactorType has to be a localizable function!");

public:
  typedef typename DiffusionFactorImp::EntityType EntityType;
  typedef typename DiffusionFactorImp::DomainFieldType DomainFieldType;
  static const size_t dimDomain = DiffusionFactorImp::dimDomain;

  // we need to distinguish three cases here (since DiffusionTensorImp may be void):
  // given two functions, a factor and a tensor
  // given a single factor (then set the tensor to default)
  // given a single tensor (then set the factor to default)
private:
  typedef EntityType E;
  typedef DomainFieldType D;
  static const size_t d = dimDomain;

  template <bool factor, bool tensor, bool anything = true>
  struct Helper
  {
    static_assert(AlwaysFalse<DiffusionTensorImp>::value, "Unsupported combination of functions given!");
  };

  // given both
  template <bool anything>
  struct Helper<false, false, anything>
  {
    static_assert(XT::Functions::is_localizable_function<DiffusionTensorImp>::value,
                  "DiffusionTensorType has to be a localizable function!");
    static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                  "EntityTypes have to agree!");
    static_assert(
        std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
        "DomainFieldTypes have to agree!");
    static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions have to agree!");
    static_assert(DiffusionFactorImp::dimRange == 1, "DiffusionFactorType has to be scalar!");
    static_assert(DiffusionFactorImp::dimRangeCols == 1, "DiffusionFactorType has to be scalar!");
    static_assert(DiffusionTensorImp::dimRange == DiffusionTensorImp::dimDomain,
                  "DiffusionTensorType has to be matrix valued!");
    static_assert(DiffusionTensorImp::dimRangeCols == DiffusionTensorImp::dimDomain,
                  "DiffusionTensorType has to be matrix valued!");
    typedef DiffusionFactorImp FactorType;
    typedef DiffusionTensorImp TensorType;
  };

  // given only one, and this is scalar
  template <bool anything>
  struct Helper<true, false, anything>
  {
    typedef DiffusionFactorImp FactorType;
    typedef XT::Functions::ConstantFunction<E, D, d, D, d, d> TensorType;
  };

  // given only one, and this is a tensor
  template <bool anything>
  struct Helper<false, true, anything>
  {
    typedef XT::Functions::ConstantFunction<E, D, d, D, 1, 1> FactorType;
    typedef DiffusionFactorImp TensorType;
  };

  static const bool single_factor_given = DiffusionFactorImp::dimRange == 1
                                          && DiffusionFactorImp::dimRangeCols == DiffusionFactorImp::dimRange
                                          && std::is_same<DiffusionTensorImp, void>::value;
  static const bool single_tensor_given = DiffusionFactorImp::dimRange != 1
                                          && DiffusionFactorImp::dimRange == DiffusionFactorImp::dimRangeCols
                                          && std::is_same<DiffusionTensorImp, void>::value;

public:
  typedef typename Helper<single_factor_given, single_tensor_given>::FactorType DiffusionFactorType;
  typedef typename Helper<single_factor_given, single_tensor_given>::TensorType DiffusionTensorType;
  typedef LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType> derived_type;
  typedef std::tuple<std::shared_ptr<typename DiffusionFactorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionTensorType::LocalfunctionType>>
      LocalfunctionTupleType;
}; // class LocalEllipticIntegrandTraits


} // namespace internal


/**
 * \brief Computes an elliptic evaluation.
 */
template <class DiffusionFactorImp, class DiffusionTensorImp>
class LocalEllipticIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalEllipticIntegrandTraits<DiffusionFactorImp,
                                                                                  DiffusionTensorImp>,
                                           2>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalEllipticIntegrandTraits<DiffusionFactorImp, DiffusionTensorImp>,
                                        2>
      BaseType;
  typedef LocalEllipticIntegrand<DiffusionFactorImp, DiffusionTensorImp> ThisType;

public:
  typedef internal::LocalEllipticIntegrandTraits<DiffusionFactorImp, DiffusionTensorImp> Traits;
  typedef typename Traits::DiffusionFactorType DiffusionFactorType;
  typedef typename Traits::DiffusionTensorType DiffusionTensorType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

private:
  typedef XT::Common::ConstStorageProvider<DiffusionFactorType> DiffusionFactorProvider;
  typedef XT::Common::ConstStorageProvider<DiffusionTensorType> DiffusionTensorProvider;
  using typename BaseType::E;
  using typename BaseType::D;
  using BaseType::d;

public:
  LocalEllipticIntegrand(const DiffusionFactorType& diff_factor,
                         const DiffusionTensorType& diff_tensor,
                         const XT::Common::Parameter& param = {})
    : diffusion_factor_(diff_factor)
    , diffusion_tensor_(diff_tensor)
    , param_(param)
  {
  }

  LocalEllipticIntegrand(const DiffusionFactorType& diff_factor, const XT::Common::Parameter& param = {})
    : diffusion_factor_(diff_factor)
    , diffusion_tensor_(new DiffusionTensorType(
          XT::Functions::internal::UnitMatrix<typename DiffusionTensorType::RangeFieldType, dimDomain>::value()))
    , param_(param)
  {
  }


  template < // This disables the ctor if dimDomain == 1, since factor and tensor are then identical and the
      typename DiffusionType, //                                                            ctors ambiguous.
      typename = typename std::enable_if<(std::is_same<DiffusionType, DiffusionTensorType>::value) && (dimDomain > 1)
                                         && sizeof(DiffusionType)>::type>
  LocalEllipticIntegrand(const DiffusionType& diffusion, const XT::Common::Parameter& param = {})
    : diffusion_factor_(new DiffusionFactorType(1.))
    , diffusion_tensor_(diffusion)
    , param_(param)
  {
  }

  /**
   * \attention Due to the nature of XT::Common::ConstStorageProvider, this copy may leave you with a dead reference,
   *            if other is destructed.
   */
  LocalEllipticIntegrand(const ThisType& other)
    : diffusion_factor_(other.diffusion_factor())
    , diffusion_tensor_(other.diffusion_tensor())
    , param_(other.param_)
  {
  }

  LocalEllipticIntegrand(ThisType&& source) = default;

  /// \name Required by LocalVolumeIntegrandInterface< ..., 2 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.access().local_function(entity),
                           diffusion_tensor_.access().local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(const LocalfunctionTupleType& local_functions_tuple,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& test_base,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatz_base) const
  {
    const auto local_diffusion_factor = std::get<0>(local_functions_tuple);
    const auto local_diffusion_tensor = std::get<1>(local_functions_tuple);
    return order(*local_diffusion_factor, *local_diffusion_tensor, test_base, ansatz_base);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& local_functions_tuple,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& test_base,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatz_base,
                const Dune::FieldVector<D, d>& localPoint,
                Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor = std::get<0>(local_functions_tuple);
    const auto local_diffusion_tensor = std::get<1>(local_functions_tuple);
    evaluate(*local_diffusion_factor, *local_diffusion_tensor, test_base, ansatz_base, localPoint, ret);
  }

  /// \}
  /// \name Actual implementations of order
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(const XT::Functions::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& local_diffusion_factor,
               const XT::Functions::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& local_diffusion_tensor,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& test_base,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatz_base) const
  {
    return local_diffusion_factor.order() + local_diffusion_tensor.order()
           + std::max(ssize_t(test_base.order()) - 1, ssize_t(0))
           + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0));
  }

  /// \}
  /// \name Actual implementations of evaluate
  /// \{

  template <class R, size_t r>
  void evaluate(const XT::Functions::LocalfunctionInterface<E, D, d, R, 1, 1>& local_diffusion_factor,
                const XT::Functions::LocalfunctionInterface<E, D, d, R, d, d>& local_diffusion_tensor,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, 1>& test_base,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, 1>& ansatz_base,
                const Dune::FieldVector<D, d>& localPoint,
                Dune::DynamicMatrix<R>& ret) const
  {
    typedef XT::Common::FieldMatrix<R, d, d> TensorType;
    ret *= 0.0;
    // evaluate local functions
    const auto diffusion_factor_value = local_diffusion_factor.evaluate(localPoint, param_);
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(localPoint, param_);
    const auto diffusion_value = diffusion_tensor_value * diffusion_factor_value;
    // evaluate bases
    const auto testGradients = test_base.jacobian(localPoint, param_);
    const auto ansatzGradients = ansatz_base.jacobian(localPoint, param_);
    // compute elliptic evaluation
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        for (size_t rr = 0; rr < r; ++rr) {
          retRow[jj] += (diffusion_value * ansatzGradients[jj][rr]) * testGradients[ii][rr];
        }
      }
    }
  } // ... evaluate(...)

  /// \}
  /// \name Access to data functions (allows the evaluation to be used as traits handling multiple szenarios).
  /// \{

  const DiffusionFactorType& diffusion_factor() const
  {
    return diffusion_factor_.access();
  }

  const DiffusionTensorType& diffusion_tensor() const
  {
    return diffusion_tensor_.access();
  }

  /// \}

private:
  const DiffusionFactorProvider diffusion_factor_;
  const DiffusionTensorProvider diffusion_tensor_;
  const XT::Common::Parameter param_;
}; // class LocalEllipticIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
