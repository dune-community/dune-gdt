// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH

#include <boost/multi_array.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "../quadrature.hh"
#include "reconstructed_function.hh"
#include "slopes.hh"
#include "internal.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, class JacobianWrapperType>
class LocalLinearReconstructionOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t axis_size = 3;
  using EntityType = typename GridLayerType::template Codim<0>::Entity;
  using IndexSetType = typename GridLayerType::IndexSet;
  using DomainType = typename BoundaryValueType::DomainType;
  using DomainFieldType = typename BoundaryValueType::DomainFieldType;
  using RangeType = typename BoundaryValueType::RangeType;
  using RangeFieldType = typename BoundaryValueType::RangeFieldType;
  using Quadrature1dType = Dune::QuadratureRule<DomainFieldType, 1>;
  using IntersectionType = typename GridLayerType::Intersection;
  using IntersectionVectorType = FieldVector<IntersectionType, 2 * dimDomain>;
  using IntersectionLocalCoordType = typename IntersectionType::Geometry::LocalCoordinate;
  using AnalyticalFluxLocalfunctionType = typename AnalyticalFluxType::LocalfunctionType;
  using StateRangeType = typename AnalyticalFluxLocalfunctionType::StateRangeType;
  using BoundaryInfoType = typename XT::Grid::BoundaryInfo<IntersectionType>;
  using ReconstructedFunctionType =
      ReconstructedLocalizableFunction<GridLayerType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;
  using VectorType = typename JacobianWrapperType::VectorType;
  using MatrixType = typename JacobianWrapperType::MatrixType;
  using StencilType = boost::multi_array<boost::optional<VectorType>, dimDomain>;
  using MultiArrayType = boost::multi_array<VectorType, dimDomain>;
  using SliceType = internal::Slice<VectorType, dimDomain>;
  using CoordsType = std::array<size_t, dimDomain>;
  using SlopeType = SlopeBase<VectorType, MatrixType>;

public:
  explicit LocalLinearReconstructionOperator(const std::vector<VectorType>& source_values,
                                             const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const SlopeType& slope,
                                             const GridLayerType& grid_layer,
                                             const XT::Common::Parameter& param,
                                             const Quadrature1dType& quadrature,
                                             ReconstructedFunctionType& reconstructed_function,
                                             XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , slope_(slope)
    , grid_layer_(grid_layer)
    , param_(param)
    , quadrature_(quadrature)
    , reconstructed_function_(reconstructed_function)
    , jacobian_wrapper_(jacobian_wrapper)
  {
    param_.set("boundary", {0.});
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    static const CoordsType stencil_sizes = []() {
      CoordsType ret;
      const auto ax_size = axis_size; // avoid linker error
      ret.fill(ax_size);
      return ret;
    }();
    thread_local StencilType stencil(stencil_sizes);
    bool valid = fill_stencil(stencil, entity);
    // In a MPI parallel run, if entity is on boundary of overlap, we do not have to reconstruct
    if (!valid)
      return;
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity)) {
      const size_t index = static_cast<size_t>(intersection.indexInInside());
      intersections[index] = intersection;
    }
    const auto entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_function_.values()[entity_index];

    // get jacobian
    bool disable_reconstruction = false;
    const auto& u_entity = source_values_[entity_index];
    auto& jac = *jacobian_wrapper_;
    if (!jac.computed() || !analytical_flux_.is_affine()) {
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      try {
        jac.get_jacobian(entity, analytical_flux_, x_in_inside_coords, u_entity, param_);
        jac.compute();
      } catch (const Dune::MathError&) {
        // Failed to compute jacobian, disable reconstruction
        disable_reconstruction = true;
      }
      if (analytical_flux_.is_affine())
        jac.jacobian() = nullptr;
    }

    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (quadrature_.size() == 1) {
        if (!disable_reconstruction) {
          // no need to reconstruct in all directions, as we are only regarding the center of the face, which will
          // always
          // have the same value assigned, independent of the slope in the other directions
          std::array<size_t, dimDomain> indices;
          indices.fill(1);
          FieldVector<VectorType, axis_size> stencil_1d, stencil_1d_char;
          FieldVector<VectorType, 2> reconstructed_values;
          for (size_t ii = 0; ii < axis_size; ++ii) { // transform to characteristic variables
            indices[dd] = ii;
            stencil_1d[ii] = *stencil(indices);
            jac.apply_inverse_eigenvectors(dd, stencil_1d[ii], stencil_1d_char[ii]);
          }
          // perform the actual reconstruction
          linear_reconstruction_1d(stencil_1d, stencil_1d_char, reconstructed_values, jac.eigenvectors(dd));

          // convert back to non-characteristic variables
          auto tmp_value = reconstructed_values[0];
          jac.apply_eigenvectors(dd, tmp_value, reconstructed_values[0]);
          tmp_value = reconstructed_values[1];
          jac.apply_eigenvectors(dd, tmp_value, reconstructed_values[1]);

          // store reconstructed values
          reconstructed_values_map.emplace(intersections[2 * dd].geometryInInside().center(), reconstructed_values[0]);
          reconstructed_values_map.emplace(intersections[2 * dd + 1].geometryInInside().center(),
                                           reconstructed_values[1]);
        } else {
          reconstructed_values_map.emplace(intersections[2 * dd].geometryInInside().center(), u_entity);
          reconstructed_values_map.emplace(intersections[2 * dd + 1].geometryInInside().center(), u_entity);
        }
      } else {
        thread_local MultiArrayType reconstructed_values(stencil_sizes);
        thread_local auto tmp_multiarray = reconstructed_values;
        tmp_multiarray.resize(stencil_sizes);

        // Transform values on stencil to characteristic variables of the current coordinate direction dd
        for (size_t ii = 0; ii < stencil.num_elements(); ++ii)
          jac.apply_inverse_eigenvectors(dd, *stencil.data()[ii], tmp_multiarray.data()[ii]);

        size_t curr_dir = dd;
        size_t last_dir = curr_dir;
        VectorType tmp_value;
        CoordsType current_sizes;
        for (size_t dir = 0; dir < dimDomain; ++dir) {
          curr_dir = (dd + dir) % dimDomain;
          // Transform to characteristic variables of the current reconstruction direction.
          if (dir > 0) {
            std::copy_n(reconstructed_values.shape(), dimDomain, current_sizes.begin());
            tmp_multiarray.resize(current_sizes);
            tmp_multiarray = reconstructed_values;
            std::for_each(
                tmp_multiarray.data(), tmp_multiarray.data() + tmp_multiarray.num_elements(), [&](VectorType& value) {
                  jac.apply_eigenvectors(last_dir, value, tmp_value);
                  jac.apply_inverse_eigenvectors(curr_dir, tmp_value, value);
                });
          } // if (dir > 0)
          // perform the actual reconstruction
          const auto& curr_quadrature = dir > 0 ? quadrature_ : get_left_right_quadrature();
          linear_reconstruction(
              curr_dir, curr_quadrature, tmp_multiarray, reconstructed_values, jac.eigenvectors(curr_dir));
          last_dir = curr_dir;
        } // dir
        // convert back to non-characteristic variables
        std::for_each(reconstructed_values.data(),
                      reconstructed_values.data() + reconstructed_values.num_elements(),
                      [&](VectorType& value) {
                        tmp_value = value;
                        jac.apply_eigenvectors(last_dir, tmp_value, value);
                      });

        // Convert coordinates on face to local entity coordinates and store reconstructed values
        internal::MultiIndexProvider<MultiArrayType> multi_indices(reconstructed_values);
        for (const auto& multi_index : multi_indices) {
          IntersectionLocalCoordType quadrature_point;
          for (size_t ii = 0; ii < dimDomain; ++ii)
            if (ii != dd)
              quadrature_point[ii < dd ? ii : ii - 1] = quadrature_[multi_index[ii]].position();
          reconstructed_values_map.emplace(
              intersections[2 * dd + multi_index[dd]].geometryInInside().global(quadrature_point),
              reconstructed_values(multi_index));
        } // multi_indices
      }
    } // dd
  } // void apply_local(...)

  // stencil is in original coordinates, stencil_char in characteristic coordinates
  // returned reconstructed values are in characteristic coordinates
  void linear_reconstruction_1d(const FieldVector<VectorType, axis_size>& stencil,
                                const FieldVector<VectorType, axis_size>& stencil_char,
                                FieldVector<VectorType, 2>& reconstructed_values,
                                const MatrixType& eigenvectors)
  {
    const auto& u_entity_char = stencil_char[1];
    const auto slope_char = slope_.get(stencil, stencil_char, eigenvectors) * 0.5;
    reconstructed_values[0] = u_entity_char - slope_char;
    reconstructed_values[1] = u_entity_char + slope_char;
  } // void linear_reconstruction_1d(...)

  // This is only used in several dimensions if we use an interface quadrature different from the midpoint quadrature
  // Todo: adapt to new slopes
  void linear_reconstruction(const size_t dir,
                             const Quadrature1dType& quadrature,
                             const MultiArrayType& stencil,
                             MultiArrayType& reconstructed_values,
                             const MatrixType& /*eigenvectors*/)
  {
    DUNE_THROW(Dune::NotImplemented, "This needs to be adapted to the new slopes");
    // resize the reconstructed_values array
    CoordsType new_shape;
    std::copy_n(stencil.shape(), dimDomain, new_shape.begin());
    new_shape[dir] = quadrature.size();
    reconstructed_values.resize(new_shape);

    // get array_views corresponding to the left, center and right values
    using IndexRangeType = typename MultiArrayType::index_range;
    using RangeArrayType = XT::Common::FieldVector<IndexRangeType, dimDomain>;
    using IndicesBuilderType = internal::IndicesBuilder<RangeArrayType, dimDomain, dimDomain - 1>;

    assert(stencil.shape()[dir] == 3);
    RangeArrayType ranges;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ranges[ii] = IndexRangeType(0, static_cast<long>(stencil.shape()[ii]));
    ranges[dir] = IndexRangeType(0);
    const auto u_left = stencil[IndicesBuilderType::build(ranges)];
    ranges[dir] = IndexRangeType(1);
    const auto u_entity = stencil[IndicesBuilderType::build(ranges)];
    ranges[dir] = IndexRangeType(2);
    const auto u_right = stencil[IndicesBuilderType::build(ranges)];

    // calculate slopes
    thread_local SliceType slope;
    std::array<size_t, dimDomain - 1> slope_sizes;
    std::copy_n(u_entity.shape(), dimDomain - 1, slope_sizes.begin());
    slope.resize(slope_sizes);
    internal::MultiIndexProvider<decltype(u_left)> multi_indices(u_left);
    // Todo: adapt to new slopes
    //    for (const auto& multi_index : multi_indices) {
    // slope(multi_index) = slope_limiter_.get(u_left(multi_index), u_entity(multi_index),
    // u_right(multi_index), eigenvectors);
    //    }

    // calculate reconstructed values
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ranges[ii] = IndexRangeType(0, static_cast<long>(new_shape[ii]));
    for (size_t ii = 0; ii < quadrature.size(); ++ii) {
      ranges[dir] = IndexRangeType(static_cast<long>(ii));
      auto reconstructed_ii = reconstructed_values[IndicesBuilderType::build(ranges)];
      for (const auto& multi_index : multi_indices) {
        reconstructed_ii(multi_index) = slope(multi_index);
        reconstructed_ii(multi_index) *= (quadrature[ii].position() - 0.5);
        reconstructed_ii(multi_index) += u_entity(multi_index);
      }
    } // ii (quadrature.size())
  } // void linear_reconstruction(...)

private:
  CoordsType center(const StencilType& stencil) const
  {
    CoordsType ret;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      assert(stencil.shape()[ii] % 2 && "Center not well-defined if one of the axis_sizes is even!");
      ret[ii] = stencil.shape()[ii] / 2;
    }
    return ret;
  }

  bool fill_stencil(StencilType& stencil, const EntityType& entity)
  {
    const int dir = -2;
    auto coords = center(stencil);
    std::fill_n(stencil.data(), stencil.num_elements(), boost::none);
    return fill_impl(stencil, entity, dir, coords);
  } // void fill(...)

private:
  bool fill_impl(StencilType& stencil, const EntityType& entity, const int dir, const CoordsType& coords)
  {
    bool ret = true;
    const auto entity_index = grid_layer_.indexSet().index(entity);
    stencil(coords) = source_values_[entity_index];
    std::vector<int> boundary_dirs;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity)) {
      const auto new_dir = intersection.indexInInside();
      if (direction_allowed(dir, new_dir) && !end_of_stencil(stencil, new_dir, coords)) {
        auto new_coords = coords;
        if (intersection.boundary() && !intersection.neighbor()) { // boundary intersections
          boundary_dirs.push_back(new_dir);
          auto boundary_value = boundary_values_.local_function(entity)->evaluate(
              intersection, entity.geometry().local(intersection.geometry().center()), source_values_[entity_index]);
          while (!end_of_stencil(stencil, new_dir, new_coords)) {
            next_coords_in_dir(new_dir, new_coords);
            stencil(new_coords) = boundary_value;
          }
        } else if (intersection.neighbor()) { // inner and periodic intersections
          const auto& outside = intersection.outside();
          next_coords_in_dir(new_dir, new_coords);
          ret = ret && fill_impl(stencil, outside, new_dir, new_coords);
        } else if (!intersection.neighbor() && !intersection.boundary()) { // processor boundary
          return false;
        }
      } // if (!end_of_stencil(...))
    } // intersections

    assert(boundary_dirs.size() <= dimDomain);
    if (boundary_dirs.size() > 1) {
      auto new_coords = coords;
      next_coords_in_dir(boundary_dirs[0], new_coords);
      const auto& boundary_value = stencil(new_coords);

      std::for_each(stencil.data(),
                    stencil.data() + stencil.num_elements(),
                    [&boundary_value](boost::optional<VectorType>& value) {
                      if (!value)
                        value = boundary_value;
                    });
    } // if (boundary_dirs.size() > 1)
    return ret;
  }

  //  get next coords in direction dir (increase or decrease coords in that direction)
  static void next_coords_in_dir(const int dir, CoordsType& coords)
  {
    assert(dir >= 0.);
    size_t udir = static_cast<size_t>(dir);
    udir % 2 ? coords[udir / 2]++ : coords[udir / 2]--;
  }

  // Direction is allowed if end of stencil is not reached and direction is not visited by another iterator.
  // Iterators never change direction, they may only spawn new iterators in the directions that have a higher
  // index (i.e. iterators that walk in x direction will spawn iterators going in y and z direction,
  // iterators going in y direction will only spawn iterators in z-direction and z iterators only walk
  // without emitting new iterators).
  static bool direction_allowed(const int dir, const int new_dir)
  {
    return new_dir == dir || new_dir / 2 > dir / 2;
  }

  static bool end_of_stencil(const StencilType& stencil, const int dir, const CoordsType& coords)
  {
    assert(dir >= 0.);
    size_t udir = static_cast<size_t>(dir);
    return coords[udir / 2] == stencil.shape()[udir / 2] - 1 || coords[udir / 2] == 0;
  }


  // quadrature rule containing left and right interface points
  static const Quadrature1dType& get_left_right_quadrature()
  {
    static const Quadrature1dType ret = left_right_quadrature();
    return ret;
  }

  static Quadrature1dType left_right_quadrature()
  {
    Quadrature1dType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  const std::vector<VectorType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SlopeType& slope_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const Quadrature1dType& quadrature_;
  ReconstructedFunctionType& reconstructed_function_;
  XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper_;
}; // class LocalLinearReconstructionOperator


template <class AnalyticalFluxImp, class BoundaryValueImp, class JacobianWrapperImp, class Traits>
class LinearReconstructionOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, class JacobianWrapperImp>
struct LinearReconstructionOperatorTraits
{
  using AnalyticalFluxType = AnalyticalFluxImp;
  using BoundaryValueType = BoundaryValueImp;
  using JacobianWrapperType = JacobianWrapperImp;
  using DomainFieldType = typename BoundaryValueType::DomainFieldType;
  using RangeFieldType = typename BoundaryValueType::DomainFieldType;
  using FieldType = DomainFieldType;
  using JacobianType = NoJacobian;
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  using ProductQuadratureType = QuadratureRule<DomainFieldType, dimDomain - 1>;
  using Quadrature1dType = Dune::QuadratureRule<DomainFieldType, 1>;
  using derived_type = LinearReconstructionOperator<AnalyticalFluxType,
                                                    BoundaryValueType,
                                                    JacobianWrapperType,
                                                    LinearReconstructionOperatorTraits>;
};


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class JacobianWrapperImp = internal::JacobianWrapper<AnalyticalFluxImp,
                                                               FieldMatrix<typename BoundaryValueImp::RangeFieldType,
                                                                           BoundaryValueImp::dimRange,
                                                                           BoundaryValueImp::dimRange>,
                                                               FieldVector<typename BoundaryValueImp::RangeFieldType,
                                                                           BoundaryValueImp::dimRange>>,
          class Traits =
              internal::LinearReconstructionOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, JacobianWrapperImp>>
class LinearReconstructionOperator : public OperatorInterface<Traits>
{
public:
  using AnalyticalFluxType = typename Traits::AnalyticalFluxType;
  using BoundaryValueType = typename Traits::BoundaryValueType;
  using JacobianWrapperType = typename Traits::JacobianWrapperType;
  using Quadrature1dType = typename Traits::Quadrature1dType;
  using ProductQuadratureType = typename Traits::ProductQuadratureType;
  using DomainFieldType = typename Traits::DomainFieldType;
  using RangeFieldType = typename Traits::RangeFieldType;
  using VectorType = typename JacobianWrapperType::VectorType;
  using MatrixType = typename JacobianWrapperType::MatrixType;
  using SlopeType = SlopeBase<VectorType, MatrixType>;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static constexpr size_t dimRange = Traits::dimRange;

  LinearReconstructionOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const SlopeType& slope = default_minmod_slope(),
                               const Quadrature1dType& quadrature_1d = default_1d_quadrature<DomainFieldType>(1))
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , slope_(slope)
    , quadrature_1d_(quadrature_1d)
    , product_quadrature_(product_quadrature_on_intersection<DomainFieldType, dimDomain>(quadrature_1d))
  {
  }

  const ProductQuadratureType& quadrature() const
  {
    return product_quadrature_;
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    static_assert(is_reconstructed_localizable_function<RangeType>::value,
                  "RangeType has to be derived from ReconstructedLocalizableFunction!");
    // evaluate cell averages
    const auto& grid_layer = source.space().grid_layer();
    const auto& index_set = grid_layer.indexSet();
    std::vector<VectorType> source_values(index_set.size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto entity_index = index_set.index(entity);
      const auto local_source = source.local_function(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator =
        LocalLinearReconstructionOperator<AnalyticalFluxType,
                                          BoundaryValueType,
                                          typename SourceType::SpaceType::GridLayerType,
                                          JacobianWrapperType>(source_values,
                                                               analytical_flux_,
                                                               boundary_values_,
                                                               slope_,
                                                               grid_layer,
                                                               param,
                                                               quadrature_1d_,
                                                               range,
                                                               jacobian_wrapper_);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(grid_layer);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  static SlopeType& default_minmod_slope()
  {
    static MinmodSlope<VectorType, MatrixType> minmod_slope_;
    return minmod_slope_;
  }

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SlopeType& slope_;
  const Quadrature1dType quadrature_1d_;
  const ProductQuadratureType product_quadrature_;
  mutable XT::Common::PerThreadValue<JacobianWrapperType> jacobian_wrapper_;
}; // class LinearReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
