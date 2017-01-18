// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
#define DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH

#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/local/discretefunction.hh>

namespace Dune {
namespace GDT {

namespace internal {

/// \todo drop and implement the is_... properly
class IsLocalOperator
{
};

/// \todo drop and implement the is_... properly
class IsLocalCouplingOperator
{
};

/// \todo drop and implement the is_... properly
class IsLocalBoundaryOperator
{
};

/// \todo drop and implement the is_... properly
class IsLocalVolumeTwoForm
{
};

/// \todo drop and implement the is_... properly
class IsLocalCouplingTwoForm
{
};

/// \todo drop and implement the is_... properly
class IsLocalBoundaryTwoForm
{
};


} // namespace internal


template <class Traits>
class LocalOperatorInterface : public XT::CRTPInterface<LocalOperatorInterface<Traits>, Traits>,
                               internal::IsLocalOperator
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   * \brief Applies the local operator.
   * \param source Should be a XT::Functions::LocalizableFunctionInterface or a ConstDiscreteFunction.
   */
  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(source, local_range));
  }
}; // class LocalOperatorInterface


template <class Traits>
class LocalCouplingOperatorInterface : public XT::CRTPInterface<LocalCouplingOperatorInterface<Traits>, Traits>,
                                       internal::IsLocalCouplingOperator
{
public:
  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor) const
  {
    this->as_imp().apply(source, intersection, local_range_entity, local_range_neighbor);
  }
}; // class LocalCouplingOperatorInterface


template <class Traits>
class LocalBoundaryOperatorInterface : public XT::CRTPInterface<LocalBoundaryOperatorInterface<Traits>, Traits>,
                                       internal::IsLocalBoundaryOperator
{
public:
  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity) const
  {
    this->as_imp().apply(source, intersection, local_range_entity);
  }
}; // class LocalBoundaryOperatorInterface


template <class Traits>
class LocalVolumeTwoFormInterface : public XT::CRTPInterface<LocalVolumeTwoFormInterface<Traits>, Traits>,
                                    internal::IsLocalVolumeTwoForm
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief Applies the local operator as a two-form.
   *  \tparam T       Traits of the test XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam A       Traits of the ansatz XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {test_base,ansatz_base}
   *  \tparam rC{T,a} dimRangeCols of the {test_base,ansatz_base}
   */
  template <class T, class A, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void apply2(const XT::Functions::LocalfunctionSetInterface<T, D, d, R, rT, rCT>& test_base,
              const XT::Functions::LocalfunctionSetInterface<A, D, d, R, rA, rCA>& ansatz_base,
              Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply2(test_base, ansatz_base, ret));
  }

  /**
   * \sa apply2
   */
  template <class T, class A, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  Dune::DynamicMatrix<R> apply2(const XT::Functions::LocalfunctionSetInterface<T, D, d, R, rT, rCT>& test_base,
                                const XT::Functions::LocalfunctionSetInterface<A, D, d, R, rA, rCA>& ansatz_base) const
  {
    Dune::DynamicMatrix<R> ret(test_base.size(), ansatz_base.size(), 0.);
    apply2(test_base, ansatz_base, ret);
    return ret;
  }
}; // class LocalVolumeTwoFormInterface


template <class Traits>
class LocalCouplingTwoFormInterface : public XT::CRTPInterface<LocalCouplingTwoFormInterface<Traits>, Traits>,
                                      internal::IsLocalCouplingTwoForm
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief Applies the local operator associated with inner faces as a two-form.
   *  \tparam TE      Traits of the entity test XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam AE      Traits of the entity ansatz XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam TN      Traits of the neighbor test XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam AN      Traits of the neighbor ansatz XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam IntersectionType
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {test_base*,ansatz_base*}
   *  \tparam rC{T,a} dimRangeCols of the {test_base*,ansatz_base*}
   */
  template <class TE,
            class AE,
            class TN,
            class AN,
            class IntersectionType,
            class D,
            size_t d,
            class R,
            size_t rT,
            size_t rCT,
            size_t rA,
            size_t rCA>
  void apply2(const XT::Functions::LocalfunctionSetInterface<TE, D, d, R, rT, rCT>& test_base_en,
              const XT::Functions::LocalfunctionSetInterface<AE, D, d, R, rA, rCA>& ansatz_base_en,
              const XT::Functions::LocalfunctionSetInterface<TN, D, d, R, rT, rCT>& test_base_ne,
              const XT::Functions::LocalfunctionSetInterface<AN, D, d, R, rA, rCA>& ansatz_base_ne,
              const IntersectionType& intersection,
              Dune::DynamicMatrix<R>& ret_en_en,
              Dune::DynamicMatrix<R>& ret_ne_ne,
              Dune::DynamicMatrix<R>& ret_en_ne,
              Dune::DynamicMatrix<R>& ret_ne_en) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply2(test_base_en,
                                              ansatz_base_en,
                                              test_base_ne,
                                              ansatz_base_ne,
                                              intersection,
                                              ret_en_en,
                                              ret_ne_ne,
                                              ret_en_ne,
                                              ret_ne_en));
  }
}; // class LocalCouplingTwoFormInterface


template <class Traits>
class LocalBoundaryTwoFormInterface : public XT::CRTPInterface<LocalBoundaryTwoFormInterface<Traits>, Traits>,
                                      internal::IsLocalBoundaryTwoForm
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief Applies the local operator associated with boundary faces as a two-form.
   *  \tparam T       Traits of the test XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam A       Traits of the ansatz XT::Functions::LocalfunctionSetInterface implementation
   *  \tparam IntersectionType
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {test_base,ansatz_base}
   *  \tparam rC{T,a} dimRangeCols of the {test_base,ansatz_base}
   */
  template <class T,
            class A,
            class IntersectionType,
            class D,
            size_t d,
            class R,
            size_t rT,
            size_t rCT,
            size_t rA,
            size_t rCA>
  void apply2(const XT::Functions::LocalfunctionSetInterface<T, D, d, R, rT, rCT>& test_base,
              const XT::Functions::LocalfunctionSetInterface<A, D, d, R, rA, rCA>& ansatz_base,
              const IntersectionType& intersection,
              Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply2(test_base, ansatz_base, intersection, ret));
  }
}; // class LocalBoundaryTwoFormInterface


/// \todo move to type_traits.hh
template <class T>
struct is_local_operator : public std::is_base_of<internal::IsLocalOperator, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_coupling_operator : public std::is_base_of<internal::IsLocalCouplingOperator, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_boundary_operator : public std::is_base_of<internal::IsLocalBoundaryOperator, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_volume_twoform : public std::is_base_of<internal::IsLocalVolumeTwoForm, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_coupling_twoform : public std::is_base_of<internal::IsLocalCouplingTwoForm, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_boundary_twoform : public std::is_base_of<internal::IsLocalBoundaryTwoForm, T>
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
