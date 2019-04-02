// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2016 - 2018)

#ifndef DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH
#define DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/integrals.hh>

#include <dune/gdt/local/finite-elements/lagrange.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, class R = double>
class FiniteVolumeGlobalBasis : public GlobalBasisInterface<GV, r, 1, R>
{
  using ThisType = FiniteVolumeGlobalBasis<GV, r, R>;
  using BaseType = GlobalBasisInterface<GV, r, 1, R>;

public:
  using BaseType::d;
  using BaseType::rC;
  using typename BaseType::D;
  using typename BaseType::E;
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalizedType;

private:
  using FiniteElementType = LocalFiniteElementInterface<D, d, R, r>;

public:
  FiniteVolumeGlobalBasis(const ThisType&) = default;
  FiniteVolumeGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  FiniteVolumeGlobalBasis(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
    , local_finite_elements_()
  {
    this->update_after_adapt();
  }

  size_t max_size() const override final
  {
    return 1;
  }

  using BaseType::localize;

  std::unique_ptr<LocalizedType> localize() const override final
  {
    return std::make_unique<LocalizedFiniteVolumeGlobalBasis>(*this);
  }

  void update_after_adapt() override final
  {
    // nothing to do
  }

private:
  class LocalizedFiniteVolumeGlobalBasis : public LocalizedGlobalFiniteElementInterface<E, r, 1, R>
  {
    using ThisType = LocalizedFiniteVolumeGlobalBasis;
    using BaseType = LocalizedGlobalFiniteElementInterface<E, r, 1, R>;

  public:
    using typename BaseType::DerivativeRangeSelector;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::DynamicDerivativeRangeType;
    using typename BaseType::DynamicRangeType;
    using typename BaseType::ElementType;
    using typename BaseType::LocalFiniteElementType;
    using typename BaseType::RangeSelector;
    using typename BaseType::RangeType;
    using typename BaseType::SingleDerivativeRangeType;

    LocalizedFiniteVolumeGlobalBasis(const FiniteVolumeGlobalBasis<GV, r, R>& self)
      : BaseType()
      , self_(self)
    {}

    LocalizedFiniteVolumeGlobalBasis(const ThisType&) = default;
    LocalizedFiniteVolumeGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    // required by XT::Functions::ElementFunctionSet

    size_t max_size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return r;
    }

  protected:
    void post_bind(const ElementType& elemnt) override final
    {
      current_local_fe_ = XT::Common::ConstStorageProvider<LocalFiniteElementInterface<D, d, R, r, 1>>(
          self_.local_finite_elements_.get(elemnt.geometry().type(), 0));
    }

  public:
    size_t size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return r;
    }

    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    using BaseType::derivatives;
    using BaseType::evaluate;
    using BaseType::jacobians;

    void evaluate(const DomainType& /*point_in_reference_element*/,
                  std::vector<RangeType>& result,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      if (result.size() < 1)
        result.resize(1);
      result[0] = 1;
    }


    void jacobians(const DomainType& /*point_in_reference_element*/,
                   std::vector<DerivativeRangeType>& result,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      if (result.size() < 1)
        result.resize(1);
      result[0] *= 0;
    }

    void derivatives(const std::array<size_t, d>& alpha,
                     const DomainType& /*point_in_reference_element*/,
                     std::vector<DerivativeRangeType>& result,
                     const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      if (result.size() < 1)
        result.resize(1);
      for (size_t jj = 0; jj < d; ++jj)
        if (alpha[jj] == 0) {
          for (size_t ii = 0; ii < r; ++ii)
            result[0][jj] = 1;
        } else {
          DUNE_THROW(Exceptions::basis_error,
                     "arbitrary derivatives are not supported!\n\n"
                         << "alpha = " << alpha);
        }
    } // ... derivatives(...)

    // These methods are default implemented in XT::Functions::ElementFunctionSetInterface and are provided here for
    // improved performance.

    void evaluate(const DomainType& /*point_in_reference_element*/,
                  std::vector<DynamicRangeType>& result,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      if (result.size() < 1)
        result.resize(1);
      RangeSelector::ensure_size(result[0]);
      for (size_t ii = 0; ii < r; ++ii)
        result[ii] = 1;
    }

    void jacobians(const DomainType& /*point_in_reference_element*/,
                   std::vector<DynamicDerivativeRangeType>& result,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      if (result.size() < 1)
        result.resize(1);
      DerivativeRangeSelector::ensure_size(result[0]);
      result[0] *= 0;
    }

    void evaluate(const DomainType& /*point_in_reference_element*/,
                  std::vector<R>& result,
                  const size_t row,
                  const size_t col = 0,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      this->assert_correct_dims(row, col, "evaluate");
      if (result.size() < 1)
        result.resize(1);
      result[0] = 1;
    }

    void jacobians(const DomainType& /*point_in_reference_element*/,
                   std::vector<SingleDerivativeRangeType>& result,
                   const size_t row,
                   const size_t col = 0,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      this->assert_correct_dims(row, col, "evaluate");
      if (result.size() < 1)
        result.resize(1);
      result[0] *= 0;
    }

    // required by LocalizedGlobalFiniteElementInterface

    const LocalFiniteElementType& finite_element() const override final
    {
      return current_local_fe_.access();
    }

    DynamicVector<R> default_data(const GeometryType& /*geometry_type*/) const override final
    {
      return DynamicVector<R>();
    }

    DynamicVector<R> backup() const override final
    {
      return DynamicVector<R>();
    }

    void restore(const ElementType& elemnt, const DynamicVector<R>& data) override final
    {
      this->bind(elemnt);
      DUNE_THROW_IF(data.size() != 0, Exceptions::finite_element_error, "data.size() = " << data.size());
    }

    using BaseType::interpolate;

    void interpolate(const std::function<RangeType(DomainType)>& element_function,
                     const int order,
                     DynamicVector<R>& dofs) const override final
    {
      if (dofs.size() != r)
        dofs.resize(r);
      dofs = DynamicVector<R>(XT::Grid::element_integral<RangeType>(this->element(), element_function, order));
      dofs /= this->element().geometry().volume();
    } // ... interpolate(...)

  private:
    const FiniteVolumeGlobalBasis<GV, r, R>& self_;
    XT::Common::ConstStorageProvider<LocalFiniteElementInterface<D, d, R, r, 1>> current_local_fe_;
  }; // class LocalizedFiniteVolumeGlobalBasis

  const GridViewType& grid_view_;
  LocalLagrangeFiniteElementFamily<D, d, R, r> local_finite_elements_;
}; // class FiniteVolumeGlobalBasis


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH
