// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_DG_DEFAULT_HH
#define DUNE_GDT_SPACES_DG_DEFAULT_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/gdt/local/finite-elements/lagrange.hh>
#include <dune/gdt/spaces/basis/default.hh>
#include <dune/gdt/spaces/mapper/discontinuous.hh>
#include <dune/gdt/spaces/dg/interface.hh>

namespace Dune {
namespace GDT {


template <class GL, int p, class R = double>
class DiscontinuousLagrangeSpace;


namespace internal {


template <class GL, int p, class R>
class DiscontinuousLagrangeSpaceTraits
{
  static_assert(XT::Grid::is_layer<GL>::value, "");
  static_assert(0 <= p, "p-adaptive case not implemented yet!");
  using G = XT::Grid::extract_grid_t<GL>;

public:
  using derived_type = DiscontinuousLagrangeSpace<GL, p, R>;
  static const constexpr int polOrder = p;
  static const constexpr size_t dimDomain = GL::dimension;
  static const constexpr size_t dimRange = 1;
  static const constexpr size_t dimRangeCols = 1;
  static const constexpr bool continuous = false;
  using GridLayerType = GL;
  using RangeFieldType = R;
  using BackendType = double;
  static const constexpr XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const constexpr Backends backend_type{Backends::gdt};
  typedef DofCommunicationChooser<GridLayerType> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;
}; // class DiscontinuousLagrangeSpaceTraits


} // namespace internal


/**
 * The following dimensions/orders/elements are tested to work:
 *
 * - 1d: orders 0, ..., 18
 * - 2d: orders 0, ..., 10 work on simplices, cubes and mixed simplices and cubes
 * - 2d: orders 11, ..., 15 also work on simplices
 * - 3d: orders 0, ..., 7 work on simplices, cubes, prisms and mixed simplices and cubes
 * - 3d: orders 8, ..., 14 also work on simplices
 * - 3d: orders 8, 9 also work on prisms
 *
 * The following dimensions/orders/elements are tested to fail:
 *
 * - 1d: orders > 18 (basis matrix fails to invert)
 * - 2d: orders > 15 on simplices (basis matrix fails to invert)
 * - 2d: orders > 10 on cubes (basis matrix fails to invert)
 * - 3d: orders > 14 on simplices(basis matrix fails to invert)
 * - 3d: orders > 7 on cubes (basis matrix fails to invert)
 * - 3d: orders > 9 on prisms (basis matrix fails to invert)
 *
 * \sa make_lagrange_local_finite_element
 */
template <class GL, int p, class R>
class DiscontinuousLagrangeSpace
    : public DgSpaceInterface<internal::DiscontinuousLagrangeSpaceTraits<GL, p, R>, GL::dimension, 1>
{
public:
  using Traits = internal::DiscontinuousLagrangeSpaceTraits<GL, p, R>;

private:
  using BaseType = DgSpaceInterface<Traits, GL::dimension, 1>;
  using ThisType = DiscontinuousLagrangeSpace<GL, p, R>;
  using D = typename GL::ctype;
  static const constexpr size_t d = BaseType::dimDomain;
  using DofCommunicationChooserType = typename Traits::DofCommunicationChooserType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::MapperType;
  using typename BaseType::FiniteElementType;
  using typename BaseType::GlobalBasisType;
  using DomainType = FieldVector<D, d>;
  using DofCommunicatorType = typename Traits::DofCommunicatorType;

private:
  using MapperImplementation = DiscontinuousMapper<GridLayerType, FiniteElementType>;
  using GlobalBasisImplementation = DefaultGlobalBasis<GridLayerType, 1, 1, R>;

public:
  DiscontinuousLagrangeSpace(GridLayerType grd_lr)
    : grid_layer_(grd_lr)
    , communicator_(DofCommunicationChooserType::create(grid_layer_))
    , backend_(0)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    // create finite elements
    for (auto&& geometry_type : grid_layer_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_lagrange_local_finite_element<D, d, R>(geometry_type, p)));
    // create mapper and basis
    mapper_ = std::make_shared<MapperImplementation>(grid_layer_, finite_elements_);
    basis_ = std::make_shared<GlobalBasisImplementation>(grid_layer_, finite_elements_);
  }

  DiscontinuousLagrangeSpace(const ThisType&) = default;
  DiscontinuousLagrangeSpace(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  const double& backend() const
  {
    return backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  const GlobalBasisType& basis() const
  {
    return *basis_;
  }

  DofCommunicatorType& dof_communicator() const
  {
    DofCommunicationChooserType::prepare(*this, *communicator_);
    return *communicator_;
  }

  const std::vector<DomainType>& lagrange_points(const EntityType& entity) const
  {
    return get_finite_element(entity.geometry().type()).lagrange_points();
  }

private:
  const FiniteElementType& get_finite_element(const GeometryType& geometry_type) const
  {
    const auto finite_element_search_result = finite_elements_->find(geometry_type);
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   entity.geometry().type() = "
                     << geometry_type);
    return *finite_element_search_result->second;
  } // ... get_finite_element(...)

  const GridLayerType grid_layer_;
  mutable std::shared_ptr<DofCommunicatorType> communicator_;
  const double backend_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<MapperImplementation> mapper_;
  std::shared_ptr<GlobalBasisImplementation> basis_;
}; // class DiscontinuousLagrangeSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DG_DEFAULT_HH
