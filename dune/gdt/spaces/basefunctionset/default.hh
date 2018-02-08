// This file is part of the dune-gdt project:
//   https://github.com/dune-comparamnity/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_DEFAULT_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_DEFAULT_HH

#include <dune/gdt/local/finite-elements/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


// forwards, required for the traits
template <class E, class R, size_t r = 1, size_t rC = 1, class F = R>
class DefaultGlobalBasis;


namespace internal {


template <class E, class R, size_t r, size_t rC, class F>
class DefaultGlobalBasisTraits
{
public:
  using derived_type = DefaultGlobalBasis<E, R, r, rC, F>;
  using EntityType = E;
  using BackendType = LocalFiniteElementInterface<typename E::Geometry::ctype, E::dimension, R, r, rC, F>;
};


} // namespace internal


template <class E, class R, size_t r, size_t rC, class F>
class DefaultGlobalBasis : public BaseFunctionSetInterface<internal::DefaultGlobalBasisTraits<E, R, r, rC, F>,
                                                           typename E::Geometry::ctype,
                                                           E::dimension,
                                                           R,
                                                           r,
                                                           rC>
{
public:
  using Traits = internal::DefaultGlobalBasisTraits<E, R, r, rC, F>;

private:
  using BaseType = BaseFunctionSetInterface<internal::DefaultGlobalBasisTraits<E, R, r, rC, F>,
                                            typename E::Geometry::ctype,
                                            E::dimension,
                                            R,
                                            r,
                                            rC>;
  using ThisType = DefaultGlobalBasis<E, R, r, rC, F>;

public:
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  DefaultGlobalBasis(const EntityType& en, const BackendType& finite_element)
    : BaseType(en)
    , finite_element_(finite_element)
  {
  }

  DefaultGlobalBasis(const ThisType&) = default;
  DefaultGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const BackendType& backend() const
  {
    return finite_element_;
  }

  size_t size() const override final
  {
    return finite_element_.basis().size();
  }

  size_t order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return finite_element_.basis().order();
  }

  void evaluate(const DomainType& xx,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    assert(this->is_a_valid_point(xx));
    ret = finite_element_.basis().evaluate(xx);
  }

  std::vector<RangeType> evaluate(const DomainType& xx, const XT::Common::Parameter& /*param*/ = {}) const
  {
    assert(this->is_a_valid_point(xx));
    return finite_element_.basis().evaluate(xx);
  }

  void jacobian(const DomainType& xx,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    assert(this->is_a_valid_point(xx));
    // evaluate jacobian of shape functions
    ret = finite_element_.basis().jacobian(xx);
    // apply transformation
    const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
    auto tmp_value = ret[0][0];
    for (size_t ii = 0; ii < finite_element_.basis().size(); ++ii) {
      J_inv_T.mv(ret[ii][0], tmp_value);
      ret[ii][0] = tmp_value;
    }
  } // ... jacobian(...)

  std::vector<JacobianRangeType> jacobian(const DomainType& xx, const XT::Common::Parameter& /*param*/ = {}) const
  {
    assert(this->is_a_valid_point(xx));
    // evaluate jacobian of shape functions
    auto ret = finite_element_.basis().jacobian(xx);
    // apply transformation
    const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
    auto tmp_value = ret[0][0];
    for (size_t ii = 0; ii < finite_element_.basis().size(); ++ii) {
      J_inv_T.mv(ret[ii][0], tmp_value);
      ret[ii][0] = tmp_value;
    }
    return ret;
  } // ... jacobian(...)

private:
  const BackendType& finite_element_;
}; // class DefaultGlobalBasis


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_DEFAULT_HH
