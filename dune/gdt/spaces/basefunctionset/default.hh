// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_DEFAULT_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_DEFAULT_HH

#include "interface.hh"

namespace Dune {
namespace GDT {


// forwards, required for the traits
template <class Fe, class E, class R = double>
class ScalarBasefunctionSet;


namespace internal {


template <class Fe, class E, class R>
class ScalarBasefunctionSetTraits
{
  using LocalFunctionTraits =
      XT::Functions::LocalfunctionSetInterface<E, typename E::Geometry::ctype, E::dimension, R, 1>;

public:
  using derived_type = ScalarBasefunctionSet<E, R>;
  using EntityType = E;
  using BackendType = Fe;
};


} // namespace internal


template <class Fe, class E, class R>
class ScalarBasefunctionSet
  : public BaseFunctionSetInterface<internal::ScalarBasefunctionSetTraits<Fe, E, R>,
                                    typename E::Geometry::ctype,
                                    E::dimension,
                                    R,
                                    1>
{
public:
  using Traits = internal::ScalarBasefunctionSetTraits<Fe, E, R>;

private:
  using BaseType = BaseFunctionSetInterface<Traits, typename E::Geometry::ctype, E::dimension, R, 1>;
  using ThisType = ScalarBasefunctionSet<Fe, E, R>;

public:
  using typename BaseType::BackendType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::JacobianRangeType;
  using typename BaseType::RangeType;

  ScalarBasefunctionSet(const EntityType& en, const BackendType& finite_element)
    : BaseType(en)
    , finite_element_(finite_element)
  {}

  ScalarBasefunctionSet(const ThisType&) = default;
  ScalarBasefunctionSet(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const BackendType& backend() const
  {
    return finite_element_;
  }

  size_t size() const override final
  {
    return finite_element_.localBasis().size();
  }

  size_t order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return finite_element_.localBasis().order();
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& xx,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    assert(this->is_a_valid_point(xx));
    finite_element_.localBasis().evaluateFunction(xx, ret);
  }

  using BaseType::jacobian;

  void jacobian(const DomainType& xx,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    assert(this->is_a_valid_point(xx));
    // evaluate jacobian of shape functions
    finite_element_.localBasis().evaluateJacobian(xx, ret);
    // apply transformation
    const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
    auto tmp_value = ret[0][0];
    for (size_t ii = 0; ii < finite_element_.localBasis().size(); ++ii) {
      J_inv_T.mv(ret[ii][0], tmp_value);
      ret[ii][0] = tmp_value;
    }
  } // ... jacobian(...)

private:
  const BackendType& finite_element_;
}; // class ScalarBasefunctionSet


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_DEFAULT_HH
