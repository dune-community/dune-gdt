// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_BASEFUNCTIONSET_DUNE_FUNCTIONS_WRAPPER_HH
#define DUNE_GDT_PLAYGROUND_SPACES_BASEFUNCTIONSET_DUNE_FUNCTIONS_WRAPPER_HH

#include <dune/common/typetraits.hh>

#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {

#if HAVE_DUNE_FUNCTIONS


// forward, to be used in the traits and to allow for specialization
template <class GL, int p, class R, size_t r, size_t rC>
class DuneFunctionsBaseFunctionSetWrapper
{
  static_assert(Dune::AlwaysFalse<GL>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GL, int p, class R, size_t r, size_t rC>
class DuneFunctionsBaseFunctionSetWrapperTraits
{
  static_assert(XT::Grid::is_view<GL>::value, "We Probably need to use TemporaryGridView from dune-xt-grid!");

public:
  typedef DuneFunctionsBaseFunctionSetWrapper<GL, p, R, r, rC> derived_type;
  typedef Functions::LagrangeDGBasis<GL, p> BackendType;
  typedef typename XT::Grid::extract_entity<GL>::type EntityType;
};


} // namespace internal


template <class GL, int p, class R>
class DuneFunctionsBaseFunctionSetWrapper<GL, p, R, 1, 1>
    : public BaseFunctionSetInterface<internal::DuneFunctionsBaseFunctionSetWrapperTraits<GL, p, R, 1, 1>,
                                      typename GL::ctype,
                                      GL::dimension,
                                      R,
                                      1,
                                      1>
{
  typedef DuneFunctionsBaseFunctionSetWrapper<GL, p, R, 1, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::DuneFunctionsBaseFunctionSetWrapperTraits<GL, p, R, 1, 1>,
                                   typename GL::ctype,
                                   GL::dimension,
                                   R,
                                   1,
                                   1>
      BaseType;

public:
  typedef internal::DuneFunctionsBaseFunctionSetWrapperTraits<GL, p, R, 1, 1> Traits;

  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

private:
  typedef typename BackendType::LocalView LocalViewType;

public:
  DuneFunctionsBaseFunctionSetWrapper(std::shared_ptr<const BackendType> bcknd, const EntityType& ent)
    : BaseType(ent)
    , backend_(bcknd)
    , local_view_(backend_->localView())
    , tmp_gradients_(local_view_.size(), JacobianRangeType(0.))
  {
    local_view_.bind(this->entity());
  }

  DuneFunctionsBaseFunctionSetWrapper(const ThisType& other) = default;
  DuneFunctionsBaseFunctionSetWrapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const
  {
    return local_view_.size();
  }

  virtual size_t order() const
  {
    return local_view_.tree().finiteElement().localBasis().order();
  }

  using BaseType::evaluate;

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const
  {
    assert(this->is_a_valid_point(xx));
    assert(ret.size() >= size());
    const auto& local_basis = local_view_.tree().finiteElement().localBasis();
    local_basis.evaluateFunction(xx, ret);
  }

  using BaseType::jacobian;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const
  {
    assert(this->is_a_valid_point(xx));
    assert(ret.size() >= size());
    const auto& local_basis = local_view_.tree().finiteElement().localBasis();
    local_basis.evaluateJacobian(xx, tmp_gradients_);
    auto jacobian_inv = this->entity().geometry().jacobianInverseTransposed(xx);
    for (size_t ii = 0; ii < size(); ++ii)
      jacobian_inv.mv(tmp_gradients_[ii][0], ret[ii]);
  } // ... jacobian(...)

private:
  const std::shared_ptr<const BackendType> backend_;
  LocalViewType local_view_;
  mutable std::vector<JacobianRangeType> tmp_gradients_;
}; // class DuneFunctionsBaseFunctionSetWrapper


#else // HAVE_DUNE_FUNCTIONS


template <class GL, int p, class R, size_t r, size_t rC = 1>
class DuneFunctionsBaseFunctionSetWrapper
{
  static_assert(Dune::AlwaysFalse<GL>::value, "You are missing dune-functions!");
};


#endif // HAVE_DUNE_FUNCTIONS

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_BASEFUNCTIONSET_DUNE_FUNCTIONS_WRAPPER_HH
