// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_PRODUCT_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_PRODUCT_HH

#include <tuple>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward
template <class... BaseFunctionSetImps>
class ProductDefault;


namespace internal {


template <typename T, typename = void>
struct getDimRange
{
  static const size_t dimRange = 1;
};

template <typename T>
struct getDimRange<T, decltype(std::declval<T>().dimRange, void())>
{
  static const size_t dimRange = T::dimRange;
};

template <class FirstType, class... Types>
struct SumDimRange
{
  static const size_t dimRange = getDimRange<FirstType>::dimRange + SumDimRange<Types...>::dimRange;
};

template <class LastType>
struct SumDimRange<LastType>
{
  static const size_t dimRange = getDimRange<LastType>::dimRange;
};

template <size_t I = 0>
struct DynamicTupleGetter
{
  template <typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), size_t>::type
  sumSize(const std::tuple<TupleArgs...>& /*tuple*/)
  {
    return 0;
  }

  template <typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), size_t>::type sumSize(const std::tuple<TupleArgs...>& tuple)
  {
    return std::get<I>(tuple).size() + DynamicTupleGetter<I + 1>::sumSize(tuple);
  }

  template <typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), void>::type order(const std::tuple<TupleArgs...>& /*tuple*/,
                                                                              size_t& /*current_order*/)
  {
  }

  template <typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), void>::type order(const std::tuple<TupleArgs...>& tuple, size_t& current_order)
  {
    if (std::get<I>(tuple).order() > current_order)
      current_order = std::get<I>(tuple).order();
    DynamicTupleGetter<I + 1>::order(tuple, current_order);
  }

  template <class DomainType, class RangeType, typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), void>::type
  evaluate(const std::tuple<TupleArgs...>& /*tuple*/,
           const DomainType& /*xx*/,
           std::vector<RangeType>& /*ret*/,
           const size_t /*first_basis_func_index*/ = 0,
           const size_t /*first_range_index*/ = 0)
  {
  }

  template <class DomainType, class RangeType, typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), void>::type evaluate(const std::tuple<TupleArgs...>& tuple,
                                                     const DomainType& xx,
                                                     std::vector<RangeType>& ret,
                                                     const size_t first_basis_func_index = 0,
                                                     const size_t first_range_index = 0)
  {
    const auto factor_ret = std::get<I>(tuple).evaluate(xx);
    const auto num_basis_funcs_in_factor = factor_ret.size();
    const auto dimRangeFactor = factor_ret[0].size();
    for (size_t basis_func = 0; basis_func < num_basis_funcs_in_factor; ++basis_func) {
      const auto& basis_func_ret = factor_ret[basis_func];
      assert(basis_func_ret.size() == dimRangeFactor);
      ret[first_basis_func_index + basis_func] = RangeType(0);
      for (size_t jj = 0; jj < dimRangeFactor; ++jj)
        ret[first_basis_func_index + basis_func][first_range_index + jj] = basis_func_ret[jj];
    }
    DynamicTupleGetter<I + 1>::evaluate(
        tuple, xx, ret, first_basis_func_index + num_basis_funcs_in_factor, first_range_index + dimRangeFactor);
  }

  template <class DomainType, class JacobianRangeType, typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), void>::type
  jacobian(const std::tuple<TupleArgs...>& /*tuple*/,
           const DomainType& /*xx*/,
           std::vector<JacobianRangeType>& /*ret*/,
           const size_t /*first_basis_func_index*/ = 0,
           const size_t /*first_range_index*/ = 0)
  {
  }

  template <class DomainType, class JacobianRangeType, typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), void>::type jacobian(const std::tuple<TupleArgs...>& tuple,
                                                     const DomainType& xx,
                                                     std::vector<JacobianRangeType>& ret,
                                                     const size_t first_basis_func_index = 0,
                                                     const size_t first_range_index = 0)
  {
    const auto factor_ret = std::get<I>(tuple).jacobian(xx);
    const auto num_basis_funcs_in_factor = factor_ret.size();
    const auto dimRangeFactor = factor_ret[0].size();
    for (size_t basis_func = 0; basis_func < num_basis_funcs_in_factor; ++basis_func) {
      const auto& basis_func_ret = factor_ret[basis_func];
      assert(basis_func_ret.size() == dimRangeFactor);
      ret[first_basis_func_index + basis_func] = JacobianRangeType(0);
      for (size_t jj = 0; jj < dimRangeFactor; ++jj)
        ret[first_basis_func_index + basis_func][first_range_index + jj] = basis_func_ret[jj];
    }
    DynamicTupleGetter<I + 1>::jacobian(
        tuple, xx, ret, first_basis_func_index + num_basis_funcs_in_factor, first_range_index + dimRangeFactor);
  }
};


template <class... BaseFunctionSetImps>
class ProductDefaultTraits
{
public:
  typedef ProductDefault<BaseFunctionSetImps...> derived_type;
  typedef double BackendType;
  typedef typename std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::EntityType EntityType;
};


} // namespace internal


template <class... BaseFunctionSetImps>
class ProductDefault
    : public BaseFunctionSetInterface<internal::ProductDefaultTraits<BaseFunctionSetImps...>,
                                      typename std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::
                                          DomainFieldType,
                                      std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::dimDomain,
                                      typename std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::
                                          RangeFieldType,
                                      internal::SumDimRange<BaseFunctionSetImps...>::dimRange,
                                      1>
{
  typedef ProductDefault<BaseFunctionSetImps...> ThisType;
  typedef BaseFunctionSetInterface<internal::ProductDefaultTraits<BaseFunctionSetImps...>,
                                   typename std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::
                                       DomainFieldType,
                                   std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::dimDomain,
                                   typename std::tuple_element<0, std::tuple<BaseFunctionSetImps...>>::type::
                                       RangeFieldType,
                                   internal::SumDimRange<BaseFunctionSetImps...>::dimRange,
                                   1>
      BaseType;

public:
  typedef internal::ProductDefaultTraits<BaseFunctionSetImps...> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  using typename BaseType::DomainType;
  using BaseType::dimRange;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  ProductDefault(const EntityType& en, BaseFunctionSetImps... basefunctionsets)
    : BaseType(en)
    , basefunctionsets_(std::move(basefunctionsets)...)
  {
  }

  ProductDefault(ThisType&& source)
    : BaseType(source)
    , basefunctionsets_(std::move(source.basefunctionsets_))
  {
  }

  ProductDefault(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return 0;
  }

  virtual size_t size() const override final
  {
    return internal::DynamicTupleGetter<0>::sumSize(basefunctionsets_);
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    size_t ret = 0;
    internal::DynamicTupleGetter<0>::order(basefunctionsets_, ret);
    return ret;
  }

  void evaluate(const DomainType& xx,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = XT::Common::Parameter()) const override final
  {
    internal::DynamicTupleGetter<0>::evaluate(basefunctionsets_, xx, ret);
  } // ... evaluate(...)

  using BaseType::evaluate;

  void jacobian(const DomainType& xx,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = XT::Common::Parameter()) const override final
  {
    internal::DynamicTupleGetter<0>::jacobian(basefunctionsets_, xx, ret);
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  std::tuple<BaseFunctionSetImps...> basefunctionsets_;
}; // class ProductDefault


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_PRODUCT_HH
