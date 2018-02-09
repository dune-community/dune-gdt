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


template <class E, class R = double>
class FiniteVolumeGlobalBasis : public GlobalBasisInterface<E, 1, 1, R>
{
  using ThisType = FiniteVolumeGlobalBasis<E, R>;
  using BaseType = GlobalBasisInterface<E, 1, 1, R>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
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

  template <class GV>
  FiniteVolumeGlobalBasis(const GridView<GV>& grid_view)
    : shape_functions_(new std::map<GeometryType, ShapeFunctionSetImplementation>())
  {
    for (auto&& geometry_type : grid_view.indexSet().types(0))
      shape_functions_->emplace(geometry_type, ShapeFunctionSetImplementation());
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

  std::unique_ptr<LocalizedBasisType> localize(const ElementType& element) const override final
  {
    return std::make_unique<LocalizedFiniteVolumeGlobalBasis>(element);
  }

private:
  class LocalizedFiniteVolumeGlobalBasis : public XT::Functions::LocalfunctionSetInterface<E, D, d, R, 1, 1>
  {
    using ThisType = LocalizedFiniteVolumeGlobalBasis;
    using BaseType = XT::Functions::LocalfunctionSetInterface<E, D, d, R, 1, 1>;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;

    LocalizedFiniteVolumeGlobalBasis(const EntityType& elemnt)
      : BaseType(elemnt)
    {
    }

    LocalizedFiniteVolumeGlobalBasis(const ThisType&) = default;
    LocalizedFiniteVolumeGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    size_t size() const override final
    {
      return 1;
    }

    size_t order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    void evaluate(const DomainType& /*xx*/,
                  std::vector<RangeType>& ret,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      ret[0] = 1;
    }

    std::vector<RangeType> evaluate(const DomainType& /*xx*/,
                                    const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return {1};
    }

    void jacobian(const DomainType& /*xx*/,
                  std::vector<JacobianRangeType>& ret,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      ret[0] *= 0;
    }

    std::vector<JacobianRangeType> jacobian(const DomainType& /*xx*/,
                                            const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return {0};
    }
  }; // class LocalizedFiniteVolumeGlobalBasis

  std::shared_ptr<std::map<GeometryType, ShapeFunctionSetImplementation>> shape_functions_;
}; // class class FiniteVolumeGlobalBasis

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_FINITE_VOLUME_HH
