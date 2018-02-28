// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2017 - 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH
#define DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH

#include <dune/localfunctions/lagrange/p0/p0localbasis.hh>

#include <dune/gdt/local/finite-elements/wrapper.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class GV, class R = double>
class FiniteVolumeGlobalBasis : public GlobalBasisInterface<GV, 1, 1, R>
{
  using ThisType = FiniteVolumeGlobalBasis<GV, R>;
  using BaseType = GlobalBasisInterface<GV, 1, 1, R>;

public:
  using typename BaseType::E;
  using typename BaseType::D;
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::ShapeFunctionsType;
  using typename BaseType::LocalizedBasisType;

private:
  using ShapeFunctionSetImplementation = LocalFiniteElementBasisWrapper<P0LocalBasis<D, R, d>, D, d, R, 1, 1>;

public:
  FiniteVolumeGlobalBasis(const ThisType&) = default;
  FiniteVolumeGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  FiniteVolumeGlobalBasis(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
    , shape_functions_(new std::map<GeometryType, ShapeFunctionSetImplementation>())
  {
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      shape_functions_->emplace(geometry_type, geometry_type);
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const ShapeFunctionsType& shape_functions(const GeometryType& geometry_type) const override final
  {
    const auto search_result = shape_functions_->find(geometry_type);
    if (search_result == shape_functions_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   geometry_type = "
                     << geometry_type);
    return search_result->second;
  }

  size_t max_size() const override final
  {
    return 1;
  }

  std::unique_ptr<LocalizedBasisType> localize() const override final
  {
    return std::make_unique<LocalizedFiniteVolumeGlobalBasis>();
  }

  std::unique_ptr<LocalizedBasisType> localize(const ElementType& element) const override final
  {
    return std::make_unique<LocalizedFiniteVolumeGlobalBasis>(element);
  }

private:
  class LocalizedFiniteVolumeGlobalBasis : public XT::Functions::LocalFunctionSetInterface<E, 1, 1, R>
  {
    using ThisType = LocalizedFiniteVolumeGlobalBasis;
    using BaseType = XT::Functions::LocalFunctionSetInterface<E, 1, 1, R>;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeSelector;
    using typename BaseType::DerivativeRangeSelector;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::SingleDerivativeRangeType;
    using typename BaseType::DynamicRangeType;
    using typename BaseType::DynamicDerivativeRangeType;

    LocalizedFiniteVolumeGlobalBasis()
      : BaseType()
    {
    }

    LocalizedFiniteVolumeGlobalBasis(const EntityType& elemnt)
      : BaseType(elemnt)
    {
    }

    LocalizedFiniteVolumeGlobalBasis(const ThisType&) = default;
    LocalizedFiniteVolumeGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    size_t max_size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 1;
    }

    size_t size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 1;
    }

    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    using BaseType::evaluate;
    using BaseType::jacobians;
    using BaseType::derivatives;

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
        if (alpha[jj] == 0)
          for (size_t ii = 0; ii < r; ++ii)
            result[0][jj] = 1;
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
  std::shared_ptr<std::map<GeometryType, ShapeFunctionSetImplementation>> shape_functions_;
}; // class class FiniteVolumeGlobalBasis

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH
