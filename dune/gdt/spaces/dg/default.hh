// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   Rene Milk       (2018)

#ifndef DUNE_GDT_SPACES_DG_DEFAULT_HH
#define DUNE_GDT_SPACES_DG_DEFAULT_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/gdt/spaces/basefunctionset/default.hh>
#include <dune/gdt/spaces/mapper/default.hh>
#include <dune/gdt/spaces/dg/interface.hh>
#include <dune/gdt/spaces/fv/default.hh>

namespace Dune {
namespace GDT {


template <class GL, int p, class R = double>
class DiscontinuousLagrangeSpace;


namespace internal {


template <class GL, int p, class R>
class DiscontinuousLagrangeSpaceTraits
{
  template <int d_ = GL::dimension, int p_ = p>
  struct dim_and_polorder_match
  {
    static const constexpr bool value = true;
  };

  template <int p_>
  struct dim_and_polorder_match<1, p_>
  {
    static const constexpr bool value = std::conditional<p_ <= 18, std::true_type, std::false_type>::type::value;
  };

  template <int p_>
  struct dim_and_polorder_match<2, p_>
  {
    // simplices: up to 15
    // cubes: up to 10
    static const constexpr bool value = std::conditional<p_ <= 10, std::true_type, std::false_type>::type::value;
  };

  template <int p_>
  struct dim_and_polorder_match<3, p_>
  {
    // simplices: up to 14
    // cubes: up to 7
    // prisms: up to 9
    static const constexpr bool value = std::conditional<p_ <= 7, std::true_type, std::false_type>::type::value;
  };

  static_assert(XT::Grid::is_layer<GL>::value, "");
  static_assert(p != 0, "This should not happen (should default to the FvSpace in this case)!");
  static_assert(0 <= p, "p-adaptive case not implemented yet!");
  static_assert(dim_and_polorder_match<>::value,
                "The LagrangeLocalFiniteElement is known to fail for these combinations!");
  using G = XT::Grid::extract_grid_t<GL>;

public:
  using derived_type = DiscontinuousLagrangeSpace<GL, p, R>;
  static const constexpr int polOrder = p;
  static const constexpr size_t dimDomain = GL::dimension;
  static const constexpr size_t dimRange = 1;
  static const constexpr size_t dimRangeCols = 1;
  static const constexpr bool continuous = false;
  using GridLayerType = GL;
  using LocalFiniteElement = LagrangeLocalFiniteElement<EquidistantPointSet, dimDomain, typename GL::ctype, R>;
  using BaseFunctionSetType = ScalarBasefunctionSet<LocalFiniteElement, XT::Grid::extract_entity_t<GL>, R>;
  using MapperType = FixedOrderScalarDiscontinuousMapper<GL, LocalFiniteElement>;
  using RangeFieldType = R;
  using BackendType = double;
  static const constexpr XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const constexpr Backends backend_type{Backends::gdt};
  typedef DofCommunicationChooser<GridLayerType> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;
}; // class DiscontinuousLagrangeSpaceTraits


} // namespace internal


/**
 * \todo Drop this specialization once the local FEs are wrapped and the LagrangeLocalFiniteElement works for order 0.
 */
template <class GL, class R>
class DiscontinuousLagrangeSpace<GL, 0, R> : public FvSpace<GL, R, 1>
{
private:
  using BaseType = FvSpace<GL, R, 1>;
  using ThisType = DiscontinuousLagrangeSpace<GL, 0, R>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::GridLayerType;

  DiscontinuousLagrangeSpace(GridLayerType grd_lr)
    : BaseType(grd_lr)
  {
  }

  std::vector<DomainType> lagrange_points(const EntityType& entity) const
  {
    return {ReferenceElements<typename GL::ctype, GL::dimension>::general(entity.type()).position(0, 0)};
  }
}; // class DiscontinuousLagrangeSpace<..., 0, ...>


/**
 * The following dimensions/orders/elements are tested to work:
 *
 * - 1d: orders 0, ..., 18
 * - 2d: orders 0, ..., 10 work on simplices, cubes and mixed simplices and cubes
 * - 2d: orders 11, ..., 15 also work on simplices but are disabled
 * - 3d: orders 0, ..., 7 work on simplices, cubes, prisms and mixed simplices and cubes
 * - 3d: orders 8, ..., 14 also work on simplices but are disabled
 * - 3d: orders 8, 9 also work on prisms but are disabled
 *
 * The following dimensions/orders/elements are tested to fail:
 *
 * - 1d: orders > 18 (basis matrix fails to invert)
 * - 2d: orders > 15 on simplices (basis matrix fails to invert)
 * - 2d: orders > 10 on cubes (basis matrix fails to invert)
 * - 3d: orders > 14 on simplices(basis matrix fails to invert)
 * - 3d: orders > 7 on cubes (basis matrix fails to invert)
 * - 3d: orders > 9 on prisms (basis matrix fails to invert)
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
  using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet, d, D, R>;
  using DofCommunicationChooserType = typename Traits::DofCommunicationChooserType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::MapperType;
  using typename BaseType::BaseFunctionSetType;
  using DomainType = typename BaseFunctionSetType::DomainType;
  using DofCommunicatorType = typename Traits::DofCommunicatorType;

  DiscontinuousLagrangeSpace(GridLayerType grd_lr)
    : grid_layer_(grd_lr)
    , communicator_(DofCommunicationChooserType::create(grid_layer_))
    , backend_(0)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , lagrange_points_(new std::map<GeometryType, std::vector<DomainType>>())
    , mapper_(nullptr)
  {
    // create finite elements and lagrange points
    for (auto&& geometry_type : grid_layer_.indexSet().types(0)) {
      //      if (geometry_type == GeometryType(GeometryType::pyramid, 3))
      //        DUNE_THROW(space_error, "Discontinuous Lagrange space does not seem to have working jacobians on pyramid
      //        grids!");
      auto fe = std::make_shared<FiniteElementType>(geometry_type, p);
      const auto& lp = fe->localInterpolation().lagrangePoints();
      std::vector<DomainType> lagrange_points(lp.size());
      for (unsigned int ii = 0; ii < lp.size(); ++ii)
        lagrange_points[ii] = lp[ii].point();
      lagrange_points_->insert(std::make_pair(geometry_type, std::move(lagrange_points)));
      finite_elements_->insert(std::make_pair(geometry_type, std::move(fe)));
    }
    // check
    //    if (d == 3 && finite_elements_->size() != 1)
    //      DUNE_THROW(space_error, "Discontinuous Lagrange space with multiple finite elements in 3d not supported
    //      (yet)!");
    // create mapper
    mapper_ = std::make_shared<MapperType>(grid_layer_, finite_elements_);
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

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    const auto finite_element_search_result = finite_elements_->find(entity.geometry().type());
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   entity.geometry().type() = "
                     << entity.geometry().type());
    const auto& finite_element = *finite_element_search_result->second;
    return BaseFunctionSetType(entity, finite_element);
  }

  DofCommunicatorType& dof_communicator() const
  {
    DofCommunicationChooserType::prepare(*this, *communicator_);
    return *communicator_;
  }

  std::vector<DomainType> lagrange_points(const EntityType& entity) const
  {
    const auto lagrange_points_search_result = lagrange_points_->find(entity.geometry().type());
    if (lagrange_points_search_result == lagrange_points_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   entity.geometry().type() = "
                     << entity.geometry().type());
    return lagrange_points_search_result->second;
  }

private:
  const GridLayerType grid_layer_;
  mutable std::shared_ptr<DofCommunicatorType> communicator_;
  const double backend_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<std::map<GeometryType, std::vector<DomainType>>> lagrange_points_;
  std::shared_ptr<MapperType> mapper_;
}; // class DiscontinuousLagrangeSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DG_DEFAULT_HH
