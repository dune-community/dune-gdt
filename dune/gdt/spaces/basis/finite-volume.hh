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
  using typename BaseType::LocalizedBasisType;
  using typename BaseType::ShapeFunctionsType;

private:
  using FiniteElementType = LocalFiniteElementInterface<D, d, R, r>;

public:
  FiniteVolumeGlobalBasis(const ThisType&) = default;
  FiniteVolumeGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  FiniteVolumeGlobalBasis(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
  {
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_local_lagrange_finite_element<D, d, R, r>(geometry_type, 0)));
  }

  const GridViewType& grid_view() const override
  {
    return grid_view_;
  }

  const ShapeFunctionsType& shape_functions(const GeometryType& geometry_type) const override final
  {
    const auto search_result = finite_elements_->find(geometry_type);
    if (search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   geometry_type = " << geometry_type);
    return search_result->second->basis();
  }

  size_t max_size() const override final
  {
    return 1;
  }

  using BaseType::localize;

  std::unique_ptr<LocalizedBasisType> localize() const override final
  {
    return std::make_unique<LocalizedFiniteVolumeGlobalBasis>();
  }

  /**
   * \note In general, we would have to check for newly created GeometryTypes and to recreate the local FEs accordingly.
   *       This is postponed until we have the LocalFiniteElementFamily.
   */
  void update_after_adapt() override final
  {
    // there is nothing to do
  }

private:
  class LocalizedFiniteVolumeGlobalBasis : public XT::Functions::ElementFunctionSetInterface<E, r, 1, R>
  {
    using ThisType = LocalizedFiniteVolumeGlobalBasis;
    using BaseType = XT::Functions::ElementFunctionSetInterface<E, r, 1, R>;

  public:
    using typename BaseType::DerivativeRangeSelector;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::DynamicDerivativeRangeType;
    using typename BaseType::DynamicRangeType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeSelector;
    using typename BaseType::RangeType;
    using typename BaseType::SingleDerivativeRangeType;

    LocalizedFiniteVolumeGlobalBasis()
      : BaseType()
    {}

    LocalizedFiniteVolumeGlobalBasis(const ThisType&) = default;
    LocalizedFiniteVolumeGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    size_t max_size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return r;
    }

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

    /**
     * \name ``These methods are required by XT::Functions::LocalizableFunctionSetInterface.''
     * \{
     **/

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

    /**
     * \}
     * \name ``These methods are default implemented in XT::Functions::LocalizableFunctionSetInterface and are
     *         overridden for improved performance.''
     * \{
     **/

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

    /**
     * \}
     * \name ``These methods (used to access an individual range dimension) are default implemented in
     *         XT::Functions::LocalizableFunctionSetInterface and are implemented for improved performance.''
     * \{
     **/

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

    /// \}

  }; // class LocalizedFiniteVolumeGlobalBasis

  const GridViewType& grid_view_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
}; // class class FiniteVolumeGlobalBasis

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH
