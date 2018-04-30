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

#include <dune/xt/grid/walker.hh>

#include <dune/xt/la/algorithms/triangular_solves.hh>
#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "../quadrature.hh"
#include "reconstructed_function.hh"
#include "slopelimiters.hh"

namespace Dune {
namespace GDT {
namespace internal {


template <class MatrixType, class VectorType, size_t dimRange, size_t num_jacobians>
class JacobianWrapper
{
protected:
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using V = XT::Common::VectorAbstraction<VectorType>;
  using EigenSolverType = typename XT::LA::EigenSolver<MatrixType>;
  using EigenSolverOptionsType = typename XT::LA::EigenSolverOptions<MatrixType>;
  using MatrixInverterOptionsType = typename XT::LA::MatrixInverterOptions<MatrixType>;

public:
  using JacobianType = typename XT::Common::FieldVector<MatrixType, num_jacobians>;

  JacobianWrapper()
    : eigenvectors_(M::create(dimRange, dimRange))
    , eigenvalues_(V::create(dimRange))
    , QR_(M::create(dimRange, dimRange))
    , tau_(V::create(dimRange))
    , jacobian_(std::make_unique<JacobianType>(eigenvectors_))
    , computed_(false)
  {
  }

  void compute(const size_t dd)
  {
    static auto eigensolver_options = create_eigensolver_options();
    const auto eigensolver = EigenSolverType((*jacobian_)[dd], eigensolver_options);
    eigenvectors_[dd] = eigensolver.real_eigenvectors();
    eigenvalues_[dd] = eigensolver.real_eigenvalues();
    QR_[dd] = eigenvectors_[dd];
    XT::LA::qr(QR_[dd], tau_[dd], permutations_[dd]);
    computed_[dd] = true;
  }

  void compute()
  {
    for (size_t dd = 0; dd < num_jacobians; ++dd)
      compute(dd);
  }

  bool computed(const size_t dd) const
  {
    return computed_[dd];
  }

  bool computed() const
  {
    for (size_t dd = 0; dd < num_jacobians; ++dd)
      if (!computed(dd))
        return false;
    return true;
  }

  std::unique_ptr<JacobianType>& jacobian()
  {
    return jacobian_;
  }

  const std::unique_ptr<JacobianType>& jacobian() const
  {
    return jacobian_;
  }

  MatrixType& jacobian(const size_t dd)
  {
    return (*jacobian_)[dd];
  }

  const MatrixType& jacobian(const size_t dd) const
  {
    return (*jacobian_)[dd];
  }

  template <class VecType>
  void apply_eigenvectors(const size_t dd, const VecType& x, VecType& ret) const
  {
    eigenvectors_[dd].mv(x, ret);
  }

  template <class VecType>
  void apply_inverse_eigenvectors(const size_t dd, const VecType& x, VecType& ret) const
  {
    VecType work = XT::Common::VectorAbstraction<VecType>::create(dimRange);
    XT::LA::solve_qr_factorized(QR_[dd], tau_[dd], permutations_[dd], ret, x, &work);
  }

protected:
  static XT::Common::Configuration create_eigensolver_options()
  {
    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[0]);
    //    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
    eigensolver_options["assert_eigendecomposition"] = "1e-6";
    eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
    eigensolver_options["disable_checks"] =
#ifdef NDEBUG
        "true";
#else
        "false";
#endif
    XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
    matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
    matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
    eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
    return eigensolver_options;
  } // ... create_eigensolver_options()

  FieldVector<MatrixType, num_jacobians> eigenvectors_;
  FieldVector<VectorType, num_jacobians> eigenvalues_;
  FieldVector<MatrixType, num_jacobians> QR_;
  FieldVector<VectorType, num_jacobians> tau_;
  FieldVector<FieldVector<int, dimRange>, num_jacobians> permutations_;
  std::unique_ptr<JacobianType> jacobian_;
  FieldVector<bool, num_jacobians> computed_;
};


} // namespace internal


template <class MultiArrayType>
class MultiIndexIterator
{
public:
  static constexpr size_t dim = MultiArrayType::dimensionality;
  static constexpr size_t array_dim = dim == 0 ? 1 : dim;
  using IndicesType = std::array<size_t, array_dim>;

  MultiIndexIterator(const MultiArrayType& multi_array, const IndicesType& indices)
    : multi_array_(multi_array)
    , indices_(indices)
  {
  }

  IndicesType& operator*()
  {
    return indices_;
  }

  const IndicesType& operator*() const
  {
    return indices_;
  }

  MultiIndexIterator& operator++()
  {
    size_t ii = 0;
    for (; ii < array_dim; ++ii) {
      indices_[ii]++;
      if (indices_[ii] < (dim == 0 ? 1 : multi_array_.shape()[ii]))
        break;
      indices_[ii] = 0;
    }
    if (ii == array_dim)
      std::copy_n(multi_array_.shape(), array_dim, indices_.begin());
    return *this;
  } // ... operator++()

  MultiIndexIterator operator++(int)
  {
    MultiIndexIterator ret = *this;
    this->operator++();
    return ret;
  } // ... operator++(int)

  bool operator==(const MultiIndexIterator& other) const
  {
    return indices_ == other.indices_;
  }

  bool operator!=(const MultiIndexIterator& other) const
  {
    return indices_ != other.indices_;
  }

private:
  const MultiArrayType& multi_array_;
  IndicesType indices_;
}; // class MultiIndexIterator<...>

template <class MultiArrayType>
class MultiIndexProvider
{
public:
  using IteratorType = MultiIndexIterator<MultiArrayType>;
  using IndicesType = typename IteratorType::IndicesType;

  MultiIndexProvider(const MultiArrayType& multi_array)
    : multi_array_(multi_array)
  {
  }

  IteratorType begin()
  {
    static const IndicesType indices = []() {
      IndicesType ret;
      ret.fill(0);
      return ret;
    }();
    return IteratorType(multi_array_, indices);
  }

  IteratorType end()
  {
    IndicesType indices{1}; // for dimension 0, the end iterator has index 1
    std::copy_n(multi_array_.shape(), MultiArrayType::dimensionality, indices.begin());
    return IteratorType(multi_array_, indices);
  }

private:
  const MultiArrayType& multi_array_;
}; // class MultiIndexProvider<...>

// Helper functor to build indices.
template <typename RangeArrayType, size_t num_ranges, size_t dimension>
struct IndicesBuilder
{
  static boost::detail::multi_array::index_gen<num_ranges, dimension> build(const RangeArrayType& ranges)
  {
    return boost::detail::multi_array::index_gen<num_ranges, dimension>(
        IndicesBuilder<RangeArrayType, num_ranges - 1, dimension>::build(ranges), ranges[num_ranges - 1]);
  }
};

// Helper functor specialization to terminate recursion.
template <typename RangeArrayType, size_t dimension>
struct IndicesBuilder<RangeArrayType, 0, dimension>
{
  static boost::detail::multi_array::index_gen<0, dimension> build(const RangeArrayType& /*ranges*/)
  {
    return boost::detail::multi_array::index_gen<0, dimension>();
  }
};

template <class RangeType, size_t dimDomain>
class Slice : public boost::multi_array<RangeType, dimDomain - 1>
{
};

template <class RangeType>
class Slice<RangeType, 1>
{
public:
  template <class MultiIndexType>
  RangeType& operator()(const MultiIndexType&)
  {
    return value_;
  }

  template <class IndicesType>
  void resize(const IndicesType&)
  {
  }

private:
  RangeType value_;
};

template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
class LocalLinearReconstructionOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t axis_size = 3;
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;
  typedef typename XT::LA::EigenSolver<MatrixType> EigenSolverType;
  typedef typename XT::LA::EigenSolverOptions<MatrixType> EigenSolverOptionsType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> Quadrature1dType;
  typedef typename GridLayerType::Intersection IntersectionType;
  typedef FieldVector<IntersectionType, 2 * dimDomain> IntersectionVectorType;
  typedef typename IntersectionType::Geometry::LocalCoordinate IntersectionLocalCoordType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  typedef typename XT::Grid::BoundaryInfo<IntersectionType> BoundaryInfoType;
  using ReconstructedFunctionType =
      ReconstructedLocalizableFunction<GridLayerType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;
  using StencilType = boost::multi_array<boost::optional<RangeType>, dimDomain>;
  using MultiArrayType = boost::multi_array<RangeType, dimDomain>;
  using SliceType = Slice<RangeType, dimDomain>;
  using CoordsType = std::array<size_t, dimDomain>;
  using JacobianWrapperType = internal::JacobianWrapper<MatrixType, std::vector<RangeFieldType>, dimRange, dimDomain>;

public:
  explicit LocalLinearReconstructionOperator(const std::vector<RangeType>& source_values,
                                             const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const GridLayerType& grid_layer,
                                             const XT::Common::Parameter& param,
                                             const Quadrature1dType& quadrature,
                                             ReconstructedFunctionType& reconstructed_function,
                                             XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(grid_layer)
    , param_(param)
    , quadrature_(quadrature)
    , reconstructed_function_(reconstructed_function)
    , jacobian_wrapper_(jacobian_wrapper)
  {
    param_.set("boundary", {0.});
  }

  void apply_local(const EntityType& entity)
  {
    auto& jac = *jacobian_wrapper_;
    static const CoordsType stencil_sizes = []() {
      CoordsType ret;
      ret.fill(axis_size);
      return ret;
    }();
    thread_local StencilType stencil(stencil_sizes);
    bool valid = fill_stencil(stencil, entity);
    // In a MPI parallel run, if entity is in overlap, we do not have to reconstruct
    if (!valid)
      return;
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity))
      intersections[intersection.indexInInside()] = intersection;
    const auto entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_function_.values()[entity_index];

    // get jacobian
    if (!jac.computed() || !analytical_flux_.is_affine()) {
      const auto& u_entity = source_values_[entity_index];
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      analytical_flux_.local_function(entity)->partial_u(x_in_inside_coords, u_entity, *jac.jacobian(), param_);
      jac.compute();
      if (analytical_flux_.is_affine())
        jac.jacobian() = nullptr;
    }

    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // Transform values on stencil to characteristic variables of the current coordinate direction dd
      thread_local MultiArrayType reconstructed_values(stencil_sizes);
      thread_local auto tmp_multiarray = reconstructed_values;
      tmp_multiarray.resize(stencil_sizes);
      for (size_t ii = 0; ii < stencil.num_elements(); ++ii)
        jac.apply_inverse_eigenvectors(dd, *stencil.data()[ii], tmp_multiarray.data()[ii]);

      size_t curr_dir = dd;
      size_t last_dir = curr_dir;
      RangeType tmp_value;
      CoordsType current_sizes;
      for (size_t dir = 0; dir < dimDomain; ++dir) {
        curr_dir = (dd + dir) % dimDomain;
        // Transform to characteristic variables of the current reconstruction direction.
        if (dir > 0) {
          std::copy_n(reconstructed_values.shape(), dimDomain, current_sizes.begin());
          tmp_multiarray.resize(current_sizes);
          tmp_multiarray = reconstructed_values;
          std::for_each(
              tmp_multiarray.data(), tmp_multiarray.data() + tmp_multiarray.num_elements(), [&](RangeType& value) {
                jac.apply_eigenvectors(last_dir, value, tmp_value);
                jac.apply_inverse_eigenvectors(curr_dir, tmp_value, value);
              });
        } // if (dir > 0)
        // perform the actual reconstruction
        const auto& curr_quadrature = dir > 0 ? quadrature_ : get_left_right_quadrature();
        linear_reconstruction(curr_dir, curr_quadrature, tmp_multiarray, reconstructed_values);
        last_dir = curr_dir;
      } // dir
      // convert back to non-characteristic variables
      std::for_each(reconstructed_values.data(),
                    reconstructed_values.data() + reconstructed_values.num_elements(),
                    [&](RangeType& value) {
                      tmp_value = value;
                      jac.apply_eigenvectors(last_dir, tmp_value, value);
                    });

      // Convert coordinates on face to local entity coordinates and store reconstructed values
      MultiIndexProvider<MultiArrayType> multi_indices(reconstructed_values);
      for (const auto& multi_index : multi_indices) {
        IntersectionLocalCoordType quadrature_point;
        for (size_t ii = 0; ii < dimDomain; ++ii)
          if (ii != dd)
            quadrature_point[ii < dd ? ii : ii - 1] = quadrature_[multi_index[ii]].position();
        reconstructed_values_map.insert(
            std::make_pair(intersections[2 * dd + multi_index[dd]].geometryInInside().global(quadrature_point),
                           reconstructed_values(multi_index)));
      } // multi_indices
    } // dd
  } // void apply_local(...)

  static void linear_reconstruction(const size_t dir,
                                    const Quadrature1dType& quadrature,
                                    const MultiArrayType& stencil,
                                    MultiArrayType& reconstructed_values)
  {
    // resize the reconstructed_values array
    CoordsType new_shape;
    std::copy_n(stencil.shape(), dimDomain, new_shape.begin());
    new_shape[dir] = quadrature.size();
    reconstructed_values.resize(new_shape);

    // get array_views corresponding to the left, center and right values
    using IndexRangeType = typename MultiArrayType::index_range;
    using RangeArrayType = XT::Common::FieldVector<IndexRangeType, dimDomain>;
    using IndicesBuilderType = IndicesBuilder<RangeArrayType, dimDomain, dimDomain - 1>;

    assert(stencil.shape()[dir] == 3);
    RangeArrayType ranges;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ranges[ii] = IndexRangeType(0, stencil.shape()[ii]);
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
    MultiIndexProvider<decltype(u_left)> multi_indices(u_left);
    for (const auto& multi_index : multi_indices) {
      const auto slope_left = u_entity(multi_index) - u_left(multi_index);
      const auto slope_center = 0.5 * (u_right(multi_index) - u_left(multi_index));
      const auto slope_right = u_right(multi_index) - u_entity(multi_index);
      slope(multi_index) = internal::ChooseLimiter<slope_limiter>::limit(slope_left, slope_right, slope_center);
    }

    // calculate reconstructed values
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ranges[ii] = IndexRangeType(0, new_shape[ii]);
    for (size_t ii = 0; ii < quadrature.size(); ++ii) {
      ranges[dir] = IndexRangeType(ii);
      auto reconstructed_ii = reconstructed_values[IndicesBuilderType::build(ranges)];
      for (const auto& multi_index : multi_indices) {
        reconstructed_ii(multi_index) = slope(multi_index);
        reconstructed_ii(multi_index) *= (quadrature[ii].position() - 0.5);
        reconstructed_ii(multi_index) += u_entity(multi_index);
      }
    } // ii (quadrature.size())
  } // static void linear_reconstruction(...)

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
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    stencil(coords) = source_values_[entity_index];
    std::vector<int> boundary_dirs;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity)) {
      const auto new_dir = intersection.indexInInside();
      if (!end_of_stencil(stencil, new_dir, coords)) {
        auto new_coords = coords;
        if (intersection.boundary() && !intersection.neighbor()) { // boundary intersections
          boundary_dirs.push_back(new_dir);
          auto boundary_value = boundary_values_.local_function(entity)->evaluate(
              intersection, entity.geometry().local(intersection.geometry().center()), source_values_[entity_index]);
          while (!end_of_stencil(stencil, new_dir, new_coords)) {
            next_coords_in_dir(new_dir, new_coords);
            stencil(new_coords) = boundary_value;
          }
        } else if (intersection.neighbor() && direction_allowed(dir, new_dir)) { // inner and periodic intersections
          const auto& outside = intersection.outside();
          next_coords_in_dir(new_dir, new_coords);
          ret = ret && fill_impl(stencil, outside, new_dir, new_coords);
        } else if (direction_allowed(dir, new_dir) && !intersection.neighbor()
                   && !intersection.boundary()) { // processor boundary
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
                    [&boundary_value](boost::optional<RangeType>& value) {
                      if (!value)
                        value = boundary_value;
                    });
    } // if (boundary_dirs.size() > 1)
    return ret;
  }

  //  get next coords in direction dir (increase or decrease coords in that direction)
  static void next_coords_in_dir(const int dir, CoordsType& coords)
  {
    dir % 2 ? coords[dir / 2]++ : coords[dir / 2]--;
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

  bool end_of_stencil(const StencilType& stencil, const int dir, const CoordsType& coords)
  {
    return coords[dir / 2] == stencil.shape()[dir / 2] - 1 || coords[dir / 2] == 0;
  }


  // quadrature rule containing left and right interface points
  static const Quadrature1dType& get_left_right_quadrature()
  {
    static Quadrature1dType ret = left_right_quadrature();
    return ret;
  }

  static Quadrature1dType left_right_quadrature()
  {
    Quadrature1dType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  const std::vector<RangeType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const Quadrature1dType quadrature_;
  ReconstructedFunctionType& reconstructed_function_;
  XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper_;
}; // class LocalLinearReconstructionOperator


#if 0
template <class SourceType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
class LocalReconstructionFvOperatorBlocked
    : public XT::Grid::Functor::Codim0<typename SourceType::SpaceType::GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  typedef typename SourceType::SpaceType SpaceType;
  typedef typename SpaceType::GridLayerType GridLayerType;
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t block_size = (dimDomain == 1) ? 2 : 4;
  static constexpr size_t num_blocks = dimRange / block_size;
  static constexpr size_t stencil_size = 2 * polOrder + 1;
  static constexpr std::array<size_t, 3> stencil = {
      {2 * polOrder + 1, dimDomain > 1 ? 2 * polOrder + 1 : 1, dimDomain > 2 ? 2 * polOrder + 1 : 1}};
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef FieldVector<RangeFieldType, block_size> LocalRangeType;
  typedef FieldMatrix<RangeFieldType, block_size, block_size> LocalMatrixType;
  typedef FieldVector<LocalMatrixType, num_blocks> MatrixType;
  typedef typename XT::LA::EigenSolver<LocalMatrixType> EigenSolverType;
  typedef typename XT::LA::EigenSolverOptions<LocalMatrixType> EigenSolverOptionsType;
  typedef typename XT::LA::MatrixInverterOptions<LocalMatrixType> MatrixInverterOptionsType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> QuadratureType;
  typedef typename GridLayerType::Intersection IntersectionType;
  typedef FieldVector<IntersectionType, 2 * dimDomain> IntersectionVectorType;
  typedef typename GridLayerType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  typedef typename Dune::FieldVector<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>
      JacobianRangeType;
  typedef typename XT::Grid::BoundaryInfo<IntersectionType> BoundaryInfoType;

public:
  explicit LocalReconstructionFvOperatorBlocked(
      SourceType& source,
      const std::vector<RangeType>& source_values,
      const AnalyticalFluxType& analytical_flux,
      const BoundaryValueType& boundary_values,
      const XT::Common::Parameter& param,
      const bool is_linear,
      const QuadratureType& quadrature,
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(source.space().grid_layer())
    , param_(param)
    , is_linear_(is_linear)
    , quadrature_(quadrature)
    , reconstructed_values_(reconstructed_values)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    param_.set("boundary", {0.});
    is_instantiated_ = true;
  }

  ~LocalReconstructionFvOperatorBlocked()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    // get cell averages on stencil
    FieldVector<int, 3> offsets(0);
    ValuesType values((FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>(
        FieldVector<RangeType, stencil[2]>(RangeType(std::numeric_limits<double>::quiet_NaN())))));
    StencilIterator::apply(source_values_, boundary_values_, values, entity, grid_layer_, -1, offsets);
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity))
      intersections[intersection.indexInInside()] = intersection;
    // get jacobians
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_values_[entity_index];
    static thread_local FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain> tau_;
    if (local_initialization_count_ != initialization_count_) {
      if (!jacobian())
        jacobian() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
      if (!eigenvectors()) {
        eigenvectors() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
        QR() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
      }
      const auto& u_entity = values[stencil[0] / 2][stencil[1] / 2][stencil[2] / 2];
      helper<dimDomain>::get_jacobian(analytical_flux_, u_entity, *(jacobian()), param_, entity);
      get_eigenvectors(*(jacobian()), *(eigenvectors()), *(QR()), tau_, permutations());
      if (is_linear_) {
        jacobian() = nullptr;
        ++local_initialization_count_;
      }
    }

    for (size_t dd = 0; dd < dimDomain; ++dd)
      helper<dimDomain>::reconstruct(dd,
                                     values,
                                     *(eigenvectors()),
                                     *(QR()),
                                     tau_,
                                     permutations(),
                                     quadrature_,
                                     reconstructed_values_map,
                                     intersections);
  } // void apply_local(...)

  static void reset()
  {
    ++initialization_count_;
  }

private:
  // quadrature rule containing left and right interface points
  static QuadratureType left_right_quadrature()
  {
    QuadratureType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper;

  static void mv(const MatrixType& Q, const RangeType& x, RangeType& ret)
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < block_size; ++ll)
        for (size_t mm = 0; mm < block_size; ++mm)
          ret[offset + ll] += Q[jj][ll][mm] * x[offset + mm];
    } // jj
  }

  static void mtv(const MatrixType& Q, const RangeType& x, RangeType& ret)
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < block_size; ++ll)
        for (size_t mm = 0; mm < block_size; ++mm)
          ret[offset + ll] += Q[jj][mm][ll] * x[offset + mm];
    } // jj
  }

  template <class anything>
  struct helper<1, anything>
  {
    static void reconstruct(size_t /*dd*/,
                            const ValuesType& values,
                            FieldVector<MatrixType, dimDomain>& eigenvectors,
                            FieldVector<MatrixType, dimDomain>& QR,
                            FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain>& tau,
                            FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations,
                            const QuadratureType& /*quadrature*/,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<RangeType, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        apply_inverse_eigenvectors(QR[0], tau[0], permutations[0], values[ii][0][0], char_values[ii]);

      // reconstruction in x direction
      FieldVector<RangeType, 2> reconstructed_values;
      // quadrature rule containing left and right interface points
      const auto left_and_right_boundary_point = get_left_right_quadrature();
      slope_reconstruction(char_values, reconstructed_values, left_and_right_boundary_point);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        // convert back to non-characteristic variables
        mv(eigenvectors[0], reconstructed_values[ii], value);
        auto quadrature_point = FieldVector<DomainFieldType, dimDomain - 1>();
        reconstructed_values_map.insert(
            std::make_pair(intersections[ii].geometryInInside().global(quadrature_point), value));
      } // ii
    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             FieldVector<MatrixType, dimDomain>& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      thread_local auto jac = std::make_unique<FieldMatrix<RangeFieldType, dimRange, dimRange>>();
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, *jac, param);
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = jj * block_size;
        for (size_t ll = 0; ll < block_size; ++ll)
          for (size_t mm = 0; mm < block_size; ++mm)
            ret[0][jj][ll][mm] = (*jac)[offset + ll][offset + mm];
      } // jj
    } // static void get_jacobian(...)
  }; // struct helper<1,...>

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             FieldVector<MatrixType, dimDomain>& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      thread_local auto jac =
          std::make_unique<FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>>();
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, *jac, param);
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        for (size_t jj = 0; jj < num_blocks; ++jj) {
          const auto offset = jj * block_size;
          for (size_t ll = 0; ll < block_size; ++ll)
            for (size_t mm = 0; mm < block_size; ++mm)
              ret[dd][jj][ll][mm] = (*jac)[dd][offset + ll][offset + mm];
        } // jj
      } // dd
    } // static void get_jacobian(...)
  }; // struct helper<3,...>

  static void
  get_eigenvectors(const FieldVector<MatrixType, dimDomain>& jacobian,
                   FieldVector<MatrixType, dimDomain>& eigenvectors,
                   FieldVector<MatrixType, dimDomain>& QR,
                   FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain>& tau,
                   FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations)
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        static XT::Common::Configuration eigensolver_options = create_eigensolver_options();
        const auto eigensolver = EigenSolverType(jacobian[dd][jj], &eigensolver_options);
        eigenvectors[dd][jj] = eigensolver.real_eigenvectors();
        QR[dd][jj] = eigenvectors[dd][jj];
        XT::LA::qr(QR[dd][jj], tau[dd][jj], permutations[dd][jj]);
      } // jj
    } // dd
  } // ... get_eigenvectors(...)

  static XT::Common::Configuration create_eigensolver_options()
  {
    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[0]);
    // XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
    eigensolver_options["assert_eigendecomposition"] = "1e-6";
    eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
    eigensolver_options["disable_checks"] = "true";
    XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
    matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
    matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
    eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
    return eigensolver_options;
  } // ... create_eigensolver_options()


  // Calculate A^{-1} x, where we have a QR decomposition with column pivoting A = QRP^T.
  // A^{-1} x = (QRP^T)^{-1} x = P R^{-1} Q^T x
  static void apply_inverse_eigenvectors(const MatrixType& QR,
                                         const FieldVector<LocalRangeType, num_blocks>& tau,
                                         const FieldVector<FieldVector<int, block_size>, num_blocks>& permutations,
                                         const RangeType& x,
                                         RangeType& ret)
  {
    LocalRangeType work;
    LocalRangeType tmp_ret, tmp_x;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      for (size_t ll = 0; ll < block_size; ++ll)
        tmp_x[ll] = x[jj * block_size + ll];
      XT::LA::solve_qr_factorized(QR[jj], tau[jj], permutations[jj], tmp_ret, tmp_x, &work);
      for (size_t ll = 0; ll < block_size; ++ll)
        ret[jj * block_size + ll] = tmp_ret[ll];
    }
  }

  const std::vector<RangeType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const bool is_linear_;
  const QuadratureType quadrature_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values_;
  //  static thread_local std::unique_ptr<JacobianRangeType> jacobian_;
  //  static thread_local std::unique_ptr<FieldVector<SparseMatrixType, dimDomain>> eigenvectors_;
  //  static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> Q_;
  //  static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> R_;
  //  static thread_local FieldVector<FieldVector<size_t, dimRange>, dimDomain> permutations_;

  // work around gcc bug 66944
  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& jacobian()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> jacobian_;
    return jacobian_;
  }

  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& eigenvectors()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> eigenvectors_;
    return eigenvectors_;
  }

  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& QR()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> QR_;
    return QR_;
  }

  static FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations()
  {
    static thread_local FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain> permutations_;
    return permutations_;
  }

  static std::atomic<size_t> initialization_count_;
  static thread_local size_t local_initialization_count_;
  static bool is_instantiated_;

}; // class LocalReconstructionFvOperatorBlocked
#endif


template <class AnalyticalFluxImp, class BoundaryValueImp, SlopeLimiters slope_limiter, class Traits>
class LinearReconstructionOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, SlopeLimiters slope_lim = SlopeLimiters::minmod>
struct LinearReconstructionOperatorTraits
{
  using AnalyticalFluxType = AnalyticalFluxImp;
  using BoundaryValueType = BoundaryValueImp;
  static const SlopeLimiters slope_limiter = slope_lim;
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
                                                    slope_limiter,
                                                    LinearReconstructionOperatorTraits>;
};


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod,
          class Traits =
              internal::LinearReconstructionOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, slope_limiter>>
class LinearReconstructionOperator : public OperatorInterface<Traits>
{
public:
  using AnalyticalFluxType = typename Traits::AnalyticalFluxType;
  using BoundaryValueType = typename Traits::BoundaryValueType;
  using Quadrature1dType = typename Traits::Quadrature1dType;
  using ProductQuadratureType = typename Traits::ProductQuadratureType;
  using DomainFieldType = typename Traits::DomainFieldType;
  using RangeFieldType = typename Traits::RangeFieldType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static constexpr size_t dimRange = Traits::dimRange;

  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using VectorType = std::vector<RangeFieldType>;
  using JacobianWrapperType = internal::JacobianWrapper<MatrixType, VectorType, dimRange, dimDomain>;

  LinearReconstructionOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const Quadrature1dType& quadrature_1d = default_1d_quadrature<DomainFieldType>(1))
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
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
    std::vector<typename SourceType::RangeType> source_values(index_set.size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto& entity_index = index_set.index(entity);
      const auto& local_source = source.local_function(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator =
        LocalLinearReconstructionOperator<AnalyticalFluxType,
                                          BoundaryValueType,
                                          typename SourceType::SpaceType::GridLayerType,
                                          slope_limiter>(source_values,
                                                         analytical_flux_,
                                                         boundary_values_,
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
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const Quadrature1dType quadrature_1d_;
  const ProductQuadratureType product_quadrature_;
  mutable XT::Common::PerThreadValue<JacobianWrapperType> jacobian_wrapper_;
}; // class LinearReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
