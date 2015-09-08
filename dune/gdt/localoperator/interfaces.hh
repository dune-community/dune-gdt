// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_INTERFACES_HH
#define DUNE_GDT_LOCALOPERATOR_INTERFACES_HH

#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/discretefunction/local.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class LocalOperatorInterface : public Stuff::CRTPInterface<LocalOperatorInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   * \brief Applies the local operator.
   * \param source Should be a Stuff::LocalizableFunctionInterface or a ConstDiscreteFunction.
   */
  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(source, local_range));
  }
}; // class LocalVolumeTwoFormInterface


template <class Traits>
class LocalVolumeTwoFormInterface : public Stuff::CRTPInterface<LocalVolumeTwoFormInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief Applies the local operator as a two-form.
   *  \tparam T       Traits of the test Stuff::LocalfunctionSetInterface implementation
   *  \tparam A       Traits of the ansatz Stuff::LocalfunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {test_base,ansatz_base}
   *  \tparam rC{T,a} dimRangeCols of the {test_base,ansatz_base}
   */
  template <class T, class A, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void apply2(const Stuff::LocalfunctionSetInterface<T, D, d, R, rT, rCT>& test_base,
              const Stuff::LocalfunctionSetInterface<A, D, d, R, rA, rCA>& ansatz_base,
              Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply2(test_base, ansatz_base, ret));
  }
}; // class LocalVolumeTwoFormInterface


namespace internal {


template <class Tt>
struct is_local_operator_helper
{
  DSC_has_typedef_initialize_once(Traits)

      static const bool is_candidate = DSC_has_typedef(Traits)<Tt>::value;
};


template <class Tt>
struct is_local_volume_twoform_helper
{
  DSC_has_typedef_initialize_once(Traits)

      static const bool is_candidate = DSC_has_typedef(Traits)<Tt>::value;
};


} // namespace internal


template <class T, bool candidate = internal::is_local_operator_helper<T>::is_candidate>
struct is_local_operator : public std::is_base_of<LocalOperatorInterface<typename T::Traits>, T>
{
};

template <class T>
struct is_local_operator<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_local_volume_twoform_helper<T>::is_candidate>
struct is_local_volume_twoform : public std::is_base_of<LocalVolumeTwoFormInterface<typename T::Traits>, T>
{
};

template <class T>
struct is_local_volume_twoform<T, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTERFACES_HH
