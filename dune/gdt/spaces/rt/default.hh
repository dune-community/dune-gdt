// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#ifndef DUNE_GDT_SPACES_RT_DEFAULT_HH
#define DUNE_GDT_SPACES_RT_DEFAULT_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/localfunctions/raviartthomas.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/gdt/spaces/basefunctionset/interface.hh>
#include <dune/gdt/spaces/mapper/default.hh>
#include <dune/gdt/spaces/rt/interface.hh>

namespace Dune {


// forwards
template <int dim, class Coordinates>
class YaspGrid;

class OneDGrid;


namespace GDT {


// forwards, required for the traits
template <class FE, class E, class R = double>
class RaviartThomasBasefunctionSet;

template <class GL, int p, class R = double>
class RaviartThomasSpace;


namespace internal {


template <class FE, class E, class R>
class RaviartThomasBasefunctionSetTraits
{
  using LocalFunctionTraits =
      XT::Functions::LocalfunctionSetInterface<E, typename E::Geometry::ctype, E::dimension, R, E::dimension>;

public:
  using derived_type = RaviartThomasBasefunctionSet<FE, E, R>;
  using EntityType = E;
  using BackendType = FE;
};


template <class GL, int p, class R>
class RaviartThomasSpaceTraits
{
  static_assert(XT::Grid::is_layer<GL>::value, "");
  static_assert(p == 0, "Not implemented yet!");
  using G = XT::Grid::extract_grid_t<GL>;

public:
  static const constexpr int polOrder = p;
  static const constexpr size_t dimDomain = GL::dimension;
  static const constexpr size_t dimRange = dimDomain;
  static const constexpr size_t dimRangeCols = 1;

private:
  template <class G_ = G, bool anything = true>
  struct ChooseLocalFiniteElement
  {
    template <class some_type_we_need_for_the_always_false = G_,
              bool single_geometry_type = Capabilities::hasSingleGeometryType<G>::v,
              bool is_simplex =
                  Capabilities::hasSingleGeometryType<G>::topologyId == Impl::SimplexTopology<dimDomain>::type::id,
              bool is_cube =
                  Capabilities::hasSingleGeometryType<G>::topologyId == Impl::CubeTopology<dimDomain>::type::id>
    struct geometry_type_switch
    {
      // This requires virtual interfaces for the local finite elements.
      static_assert(AlwaysFalse<some_type_we_need_for_the_always_false>::value,
                    "RaviartThomasSpace is only available for purely simplicial or purely cubic grids atm!");
    };

    template <class some_type_we_need_for_the_always_false>
    struct geometry_type_switch<some_type_we_need_for_the_always_false, true, true, false>
    {
      using type = RaviartThomasSimplexLocalFiniteElement<dimDomain, typename GL::ctype, R>;

      static std::shared_ptr<type> create(const Dune::GeometryType& geometry_type, const int& polorder)
      {
        return std::make_shared<type>(geometry_type, polorder);
      }
    };

    template <class some_type_we_need_for_the_always_false>
    struct geometry_type_switch<some_type_we_need_for_the_always_false, true, false, true>
    {
      using type = RaviartThomasCubeLocalFiniteElement<typename GL::ctype, R, dimDomain, p>;

      static std::shared_ptr<type> create(const Dune::GeometryType& /*geometry_type*/, const int& /*polorder*/)
      {
        return std::make_shared<type>();
      }
    };

    using type = typename geometry_type_switch<>::type;

    static std::shared_ptr<type> create(const Dune::GeometryType& geometry_type, const int& polorder)
    {
      return geometry_type_switch<>::create(geometry_type, polorder);
    }
  };

  // We need an exception for 1d YaspGrid: it reports to contain cubes but everything is a simplex in 1d.
  template <bool anything, class Coordinates>
  struct ChooseLocalFiniteElement<YaspGrid<1, Coordinates>, anything>
  {
    using type = RaviartThomasSimplexLocalFiniteElement<dimDomain, typename GL::ctype, R>;

    static std::shared_ptr<type> create(const Dune::GeometryType& geometry_type, const int& polorder)
    {
      return std::make_shared<type>(geometry_type, polorder);
    }
  };

  // We need an exception for OneDGrid: it reports to contain cubes but everything is a simplex in 1d.
  template <bool anything>
  struct ChooseLocalFiniteElement<OneDGrid, anything>
  {
    using type = RaviartThomasSimplexLocalFiniteElement<dimDomain, typename GL::ctype, R>;

    static std::shared_ptr<type> create(const Dune::GeometryType& geometry_type, const int& polorder)
    {
      return std::make_shared<type>(geometry_type, polorder);
    }
  };

public:
  using derived_type = RaviartThomasSpace<GL, p, R>;
  static const constexpr bool continuous = false;
  using GridLayerType = GL;
  using BackendType = typename ChooseLocalFiniteElement<>::type;
  using BaseFunctionSetType = RaviartThomasBasefunctionSet<BackendType, XT::Grid::extract_entity_t<GL>, R>;
  using MapperType = FixedOrderMultipleCodimMultipleGeomTypeMapper<GL, BackendType>;
  using RangeFieldType = R;
  static const constexpr XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const constexpr Backends backend_type{Backends::gdt};
  typedef DofCommunicationChooser<GridLayerType, false> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;

private:
  friend class RaviartThomasSpace<GL, p, R>;
}; // class RaviartThomasSpaceTraits


} // namespace internal


template <class FE, class E, class R>
class RaviartThomasBasefunctionSet
    : public BaseFunctionSetInterface<internal::RaviartThomasBasefunctionSetTraits<FE, E, R>,
                                      typename E::Geometry::ctype,
                                      E::dimension,
                                      R,
                                      E::dimension,
                                      1>
{
public:
  using Traits = internal::RaviartThomasBasefunctionSetTraits<FE, E, R>;

private:
  using BaseType = BaseFunctionSetInterface<Traits, typename E::Geometry::ctype, E::dimension, R, E::dimension, 1>;
  using ThisType = RaviartThomasBasefunctionSet<FE, E, R>;

public:
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using BaseType::d;

  RaviartThomasBasefunctionSet(const EntityType& en, const BackendType& finite_element, const std::vector<R>& switches)
    : BaseType(en)
    , finite_element_(finite_element)
    , switches_(switches)
  {
    if (switches_.size() != finite_element_.size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "finite_element.size() = " << finite_element_.size() << "\n   switches.size() = " << switches_.size());
  }

  RaviartThomasBasefunctionSet(const ThisType&) = default;
  RaviartThomasBasefunctionSet(ThisType&&) = default;

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

  size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return finite_element_.localBasis().order();
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& xx,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    assert(this->is_a_valid_point(xx));
    // evaluate shape functions
    finite_element_.localBasis().evaluateFunction(xx, ret);
    // flip and scale shape functions to ensure
    // - continuity of normal component and
    // - to ensure basis*integrationElementNormal = 1
    for (size_t ii = 0; ii < finite_element_.localBasis().size(); ++ii)
      ret[ii] *= switches_[ii];
    // apply piola transformation
    const auto J_T = this->entity().geometry().jacobianTransposed(xx);
    const auto det_J_T = std::abs(J_T.determinant());
    RangeType tmp_value;
    for (size_t ii = 0; ii < finite_element_.localBasis().size(); ++ii) {
      J_T.mtv(ret[ii], tmp_value);
      tmp_value /= det_J_T;
      ret[ii] = tmp_value;
    }
  } // ... evaluate(...)

  using BaseType::jacobian;

  void jacobian(const DomainType& xx,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    assert(this->is_a_valid_point(xx));
    // evaluate jacobian of shape functions
    finite_element_.localBasis().evaluateJacobian(xx, ret);
    // flip and scale shape functions to ensure
    // - continuity of normal component and
    // - to ensure basis*integrationElementNormal = 1
    for (size_t ii = 0; ii < finite_element_.localBasis().size(); ++ii)
      ret[ii] *= switches_[ii];
    // apply piola transformation
    const auto J_T = this->entity().geometry().jacobianTransposed(xx);
    const auto det_J_T = std::abs(J_T.determinant());
    const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
    auto tmp_jacobian_row = ret[0][0];
    for (size_t ii = 0; ii < finite_element_.localBasis().size(); ++ii) {
      for (size_t jj = 0; jj < d; ++jj) {
        J_inv_T.mtv(ret[ii][jj], tmp_jacobian_row);
        J_T.mtv(tmp_jacobian_row, ret[ii][jj]);
        ret[ii][jj] /= det_J_T;
      }
    }
  } // ... jacobian(...)

private:
  const BackendType& finite_element_;
  const std::vector<R>& switches_;
}; // class RaviartThomasBasefunctionSet


template <class GL, int p, class R>
class RaviartThomasSpace
    : public RtSpaceInterface<internal::RaviartThomasSpaceTraits<GL, p, R>, GL::dimension, GL::dimension>
{
public:
  using Traits = internal::RaviartThomasSpaceTraits<GL, p, R>;

private:
  using BaseType = RtSpaceInterface<Traits, GL::dimension, GL::dimension>;
  using ThisType = RaviartThomasSpace<GL, p, R>;
  using D = typename GL::ctype;
  static const constexpr size_t d = BaseType::dimDomain;
  using FiniteElementType = typename BaseType::BackendType;
  typedef typename Traits::DofCommunicationChooserType DofCommunicationChooserType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::MapperType;
  using typename BaseType::BaseFunctionSetType;
  typedef typename Traits::DofCommunicatorType DofCommunicatorType;

  RaviartThomasSpace(GridLayerType grd_lr)
    : grid_layer_(grd_lr)
    , communicator_(DofCommunicationChooserType::create(grid_layer_))
    , backend_(0)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , geometry_to_local_DoF_indices_map_(new std::map<GeometryType, std::vector<size_t>>())
    , entity_indices_(new ZeroOrderScalarDiscontinuousMapper<GL>(grid_layer_)) // <- We need unique indices for codim 0
    , switches_(new std::vector<std::vector<R>>(entity_indices_->size())) //       entities: this cannot be achieved by
    , mapper_(nullptr) //                                                     the grid layers index set for mixed grids.
    , communicator_prepared_(false)
  {
    // create finite elements
    for (auto&& geometry_type : grid_layer_.indexSet().types(0)) {
      auto finite_element = Traits::template ChooseLocalFiniteElement<>::create(geometry_type, p);
      finite_elements_->insert(std::make_pair(geometry_type, finite_element));
      // this is only valid for p0 elements
      geometry_to_local_DoF_indices_map_->insert(
          std::make_pair(geometry_type, std::vector<size_t>(finite_element->size())));
    }
    // compute local-key-to-intersection relationship
    for (const auto& geometry_type_and_finite_element_ptr : *finite_elements_) {
      const auto& geometry_type = geometry_type_and_finite_element_ptr.first;
      const auto& finite_element = *geometry_type_and_finite_element_ptr.second;
      const auto& coeffs = finite_element.localCoefficients();
      auto& local_key_to_intersection_map = geometry_to_local_DoF_indices_map_->at(geometry_type);
      for (unsigned int ii = 0; ii < coeffs.size(); ++ii) {
        const auto& local_key = coeffs.localKey(ii);
        if (local_key.index() != 0)
          DUNE_THROW(XT::Common::Exceptions::internal_error, "This must not happen for p0!");
        if (local_key.codim() != 1)
          DUNE_THROW(XT::Common::Exceptions::internal_error, "This must not happen for p0!");
        local_key_to_intersection_map[ii] = local_key.subEntity();
      }
    }
    // compute scaling to ensure basis*integrationElementNormal = 1
    std::map<GeometryType, std::vector<R>> geometry_to_scaling_factors_map;
    for (const auto& geometry_type_and_finite_element_ptr : *finite_elements_) {
      const auto& geometry_type = geometry_type_and_finite_element_ptr.first;
      const auto& finite_element = *geometry_type_and_finite_element_ptr.second;
      const auto& shape_functions = finite_element.localBasis();
      std::vector<typename BaseFunctionSetType::RangeType> shape_functions_evaluations(shape_functions.size());
      const auto& reference_element = ReferenceElements<D, d>::general(geometry_type);
      const auto num_intersections = reference_element.size(1);
      geometry_to_scaling_factors_map.insert(std::make_pair(geometry_type, std::vector<R>(num_intersections, R(1.))));
      auto& scaling_factors = geometry_to_scaling_factors_map.at(geometry_type);
      const auto& DoF_indices = geometry_to_local_DoF_indices_map_->at(geometry_type);
      for (auto&& intersection_index : XT::Common::value_range(num_intersections)) {
        const auto& normal = reference_element.integrationOuterNormal(intersection_index);
        const auto& intersection_center = reference_element.position(intersection_index, 1);
        const auto DoF_index = DoF_indices[intersection_index];
        shape_functions.evaluateFunction(intersection_center, shape_functions_evaluations);
        scaling_factors[intersection_index] /= shape_functions_evaluations[DoF_index] * normal;
      }
    }
    // compute switches (as signs of the scaling factor) to ensure continuity of the normal component
    for (auto&& entity : elements(grid_layer_)) {
      const auto geometry_type = entity.geometry().type();
      const auto& finite_element = *finite_elements_->at(geometry_type);
      const auto& coeffs = finite_element.localCoefficients();
      const auto entity_index = entity_indices_->mapToGlobal(entity, 0);
      (*switches_)[entity_index] = geometry_to_scaling_factors_map.at(geometry_type);
      auto& local_switches = switches_->at(entity_index);
      for (auto&& intersection : intersections(grid_layer_, entity)) {
        if (intersection.neighbor() && entity_index < entity_indices_->mapToGlobal(intersection.outside(), 0)) {
          const auto intersection_index = XT::Common::numeric_cast<unsigned int>(intersection.indexInInside());
          for (unsigned int ii = 0; ii < coeffs.size(); ++ii) {
            const auto& local_key = coeffs.localKey(ii);
            const auto DoF_subentity_index = local_key.subEntity();
            if (local_key.codim() == 1 && DoF_subentity_index == intersection_index)
              local_switches[DoF_subentity_index] *= -1.;
          }
        }
      }
    }
    // create mapper
    mapper_ = std::make_shared<MapperType>(grid_layer_, finite_elements_);
  } // RaviartThomasSpace(...)

  RaviartThomasSpace(const ThisType&) = default;
  RaviartThomasSpace(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  const BackendType& backend() const
  {
    // this only works because we know that there can only be one geometry.
    return *finite_elements_->begin()->second;
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
                 "This must not happen after the checks in the ctor, the grid layer did not report all geometry types!"
                 "\n   entity.geometry().type() = "
                     << entity.geometry().type());
    const auto& finite_element = *finite_element_search_result->second;
    return BaseFunctionSetType(entity, finite_element, (*switches_)[entity_indices_->mapToGlobal(entity, 0)]);
  }

  DofCommunicatorType& dof_communicator() const
  {
    if (!communicator_prepared_) {
      communicator_prepared_ = DofCommunicationChooserType::prepare(*this, *communicator_);
    }
    return *communicator_;
  }

  // this makes sense only for for p0 (and only on simplices?)
  // once we export the local finite element, this should become obsolete!
  std::vector<size_t> local_DoF_indices(const EntityType& entity) const
  {
    const auto search_result = geometry_to_local_DoF_indices_map_->find(entity.geometry().type());
    if (search_result == geometry_to_local_DoF_indices_map_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen after the checks in the ctor, the grid layer did not report all geometry types!"
                 "\n   entity.geometry().type() = "
                     << entity.geometry().type());
    return search_result->second;
  } // ... local_DoF_indices(...)

private:
  const GridLayerType grid_layer_;
  mutable std::shared_ptr<DofCommunicatorType> communicator_;
  const double backend_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<std::map<GeometryType, std::vector<size_t>>> geometry_to_local_DoF_indices_map_;
  std::shared_ptr<ZeroOrderScalarDiscontinuousMapper<GL>> entity_indices_;
  std::shared_ptr<std::vector<std::vector<R>>> switches_;
  std::shared_ptr<MapperType> mapper_;
  mutable bool communicator_prepared_;
}; // class RaviartThomasSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RT_DEFAULT_HH
