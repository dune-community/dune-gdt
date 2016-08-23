// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_OPERATORS_INTERFACES_HH
#define DUNE_GDT_OPERATORS_INTERFACES_HH

#include <type_traits>

#include <dune/common/deprecated.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/xt/la/container/interfaces.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/common/configuration.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


/**
 * \note The methods apply() and apply2() do not support threaded assembly atm, else we would have to add a use_tbb
 *       switch (to a lot of methods). Either this gets merged with the new Parameter or we need a different paradigm.
 */
template <class Traits>
class OperatorInterface : public XT::CRTPInterface<OperatorInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::FieldType FieldType;
  //  typedef typename Traits::JacobianType JacobianType;

  /// \name Methods that have to be implemented by any derived class
  /// \{

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range) const
  {
    CHECK_CRTP(this->as_imp().apply(source, range));
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& range, const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp().apply2(range, source));
    return this->as_imp().apply2(range, source);
  }

  //  template< class SourceType >
  //  JacobianType jacobian(const SourceType& source) const
  //  {
  //    CHECK_CRTP(this->as_imp().jacobian(source));
  //    return this->as_imp().jacobian(source);
  //  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& range, SourceType& source, const XT::Common::Configuration& opts) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply_inverse(range, source, opts));
  }

  std::vector<std::string> invert_options() const
  {
    CHECK_CRTP(this->as_imp().invert_options());
    return this->as_imp().invert_options();
  }

  XT::Common::Configuration invert_options(const std::string& type) const
  {
    CHECK_CRTP(this->as_imp().invert_options(type));
    return this->as_imp().invert_options(type);
  }

  /// \}
  /// \name Provided by the interface for convenience
  /// \{

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& range, SourceType& source, const std::string& type) const
  {
    apply_inverse(range, source, invert_options(type));
  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& range, SourceType& source) const
  {
    auto type = invert_options();
    apply_inverse(range, source, type.size() > 0 ? type[0] : "");
  }

  template <class RangeType>
  FieldType induced_norm(const RangeType& range) const
  {
    return std::sqrt(apply2(range, range));
  }

  /// \}
}; // class OperatorInterface


namespace internal {


template <class Tt>
struct is_operator_helper
{
  DSC_has_typedef_initialize_once(Traits)

      static const bool is_candidate = DSC_has_typedef(Traits)<Tt>::value;
};


} // namespace internal


template <class T, bool candidate = internal::is_operator_helper<T>::is_candidate>
struct is_operator : public std::is_base_of<OperatorInterface<typename T::Traits>, T>
{
};

template <class T>
struct is_operator<T, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
