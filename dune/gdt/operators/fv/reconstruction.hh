// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH

#include <dune/common/fvector.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/container/eigen.hh>

#include "eigensolver.hh"
#include "slopelimiters.hh"

namespace Dune {
namespace GDT {


template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter,
          class EigenSolverType =
              typename Dune::GDT::DefaultEigenSolver<typename AnalyticalFluxType::LocalfunctionType>>
class LocalReconstructionFvOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t stencil_size = 2 * polOrder + 1;
  static constexpr std::array<size_t, 3> stencil = {
      {2 * polOrder + 1, dimDomain > 1 ? 2 * polOrder + 1 : 1, dimDomain > 2 ? 2 * polOrder + 1 : 1}};
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef typename EigenSolverType::VectorType VectorType;
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> QuadratureType;
  typedef FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> IntersectionVectorType;
  typedef FieldVector<FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>, stencil[0]> ValuesType;
  typedef typename GridLayerType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;

public:
  explicit LocalReconstructionFvOperator(
      const std::vector<RangeType> source_values,
      const AnalyticalFluxType& analytical_flux,
      const BoundaryValueType& boundary_values,
      const GridLayerType& grid_layer,
      const XT::Common::Parameter& param,
      const bool is_linear,
      const QuadratureType& quadrature,
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(grid_layer)
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

  ~LocalReconstructionFvOperator()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    // get cell averages on stencil
    FieldVector<size_t, 3> offsets(0);
    ValuesType values;
    StencilIterator::apply(source_values_, boundary_values_, values, entity, grid_layer_, -1, offsets);
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity)) {
      const auto& n = intersection.unitOuterNormal(intersection.geometry().local(intersection.geometry().center()));
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        if (XT::Common::FloatCmp::eq(n[dd], -1.))
          intersections[2 * dd] = intersection;
        else if (XT::Common::FloatCmp::eq(n[dd], 1.))
          intersections[2 * dd + 1] = intersection;
      }
    }
    // get jacobians
    const auto& u_entity = values[stencil[0] / 2][stencil[1] / 2][stencil[2] / 2];
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_values_[entity_index];
    const auto flux_local_func = analytical_flux_.local_function(entity);
    if (!is_linear_ || !eigensolver_)
      eigensolver_ = XT::Common::make_unique<EigenSolverType>(
          *flux_local_func, entity.geometry().local(entity.geometry().center()), u_entity, param_);
    const auto& eigenvectors = eigensolver_->eigenvectors();
    const auto& eigenvectors_inverse = eigensolver_->eigenvectors_inverse();
    for (size_t dd = 0; dd < dimDomain; ++dd)
      helper<>::reconstruct(
          dd, values, eigenvectors[dd], eigenvectors_inverse[dd], quadrature_, reconstructed_values_map, intersections);
  } // void apply_local(...)

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

  template <class anything>
  struct helper<1, anything>
  {
    static void reconstruct(size_t /*dd*/,
                            const ValuesType& values,
                            const MatrixType& eigenvectors,
                            const MatrixType& eigenvectors_inverse,
                            const QuadratureType& /*quadrature*/,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<RangeType, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        eigenvectors_inverse.mv(values[ii][0][0], char_values[ii]);

      // reconstruction in x direction
      FieldVector<RangeType, 2> reconstructed_values;
      // quadrature rule containing left and right interface points
      const auto left_and_right_boundary_point = left_right_quadrature();
      slope_reconstruction(char_values, reconstructed_values, left_and_right_boundary_point);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        // convert back to non-characteristic variables
        eigenvectors.mv(reconstructed_values[ii], value);
        auto quadrature_point = FieldVector<DomainFieldType, dimDomain - 1>();
        reconstructed_values_map.insert(
            std::make_pair(intersections[ii].geometryInInside().global(quadrature_point), value));
      } // ii
    } // static void reconstruct()
  }; // struct helper<1,...>

  template <class anything>
  struct helper<2, anything>
  {

    static void reconstruct(size_t dd,
                            const ValuesType& values,
                            const MatrixType& eigenvectors,
                            const MatrixType& eigenvectors_inverse,
                            const QuadratureType& quadrature,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y-direction accordingly. For that purpose, define new
      // coordinates x', y'. First convert to characteristic variables and reorder x, y to x', y'.
      // Reordering is done such that the indices in char_values are in the order y', x'.
      FieldVector<FieldVector<RangeType, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          if (dd == 0)
            eigenvectors_inverse.mv(values[ii][jj][0], char_values[jj][ii]);
          else if (dd == 1)
            eigenvectors_inverse.mv(values[ii][jj][0], char_values[ii][jj]);
        }
      }

      // reconstruction in x' direction
      // first index: dimension 2 for left/right interface
      // second index: y' direction
      FieldVector<FieldVector<RangeType, stencil_size>, 2> x_reconstructed_values;
      std::vector<RangeType> result(2);
      const auto left_and_right_boundary_point = left_right_quadrature();
      for (size_t jj = 0; jj < stencil_size; ++jj) {
        slope_reconstruction(char_values[jj], result, left_and_right_boundary_point);
        x_reconstructed_values[0][jj] = result[0];
        x_reconstructed_values[1][jj] = result[1];
      }

      const auto& num_quad_points = quadrature.size();
      // reconstruction in y' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      FieldVector<std::vector<RangeType>, 2> reconstructed_values((std::vector<RangeType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        slope_reconstruction(x_reconstructed_values[ii], reconstructed_values[ii], quadrature);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          // convert back to non-characteristic variables and to FieldVector instead of EigenVector
          eigenvectors.mv(reconstructed_values[ii][jj], value);
          auto quadrature_point = quadrature[jj].position();
          reconstructed_values_map.insert(
              std::make_pair(intersections[2 * dd + ii].geometryInInside().global(quadrature_point), value));
        } // jj
      } // ii
    } // static void reconstruct()
  }; // helper<2,...>

  template <class anything>
  struct helper<3, anything>
  {

    static void reconstruct(size_t dd,
                            const ValuesType& values,
                            const MatrixType& eigenvectors,
                            const MatrixType& eigenvectors_inverse,
                            const QuadratureType& quadrature,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y- and z-direction accordingly. For that purpose, define new
      // coordinates x', y', z'. First convert to characteristic variables and reorder x, y, z to x', y', z'.
      // Reordering is done such that the indices in char_values are in the order z', y', x'
      FieldVector<FieldVector<FieldVector<RangeType, stencil_size>, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            if (dd == 0)
              eigenvectors_inverse.mv(values[ii][jj][kk], char_values[kk][jj][ii]);
            else if (dd == 1)
              eigenvectors_inverse.mv(values[ii][jj][kk], char_values[ii][kk][jj]);
            else if (dd == 2)
              eigenvectors_inverse.mv(values[ii][jj][kk], char_values[jj][ii][kk]);
          }
        }
      }

      // reconstruction in x' direction
      // first index: z' direction
      // second index: dimension 2 for left/right interface
      // third index: y' direction
      FieldVector<FieldVector<FieldVector<RangeType, stencil_size>, 2>, stencil_size> x_reconstructed_values;
      std::vector<RangeType> result(2);
      const auto left_and_right_boundary_point = left_right_quadrature();
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          slope_reconstruction(char_values[kk][jj], result, left_and_right_boundary_point);
          x_reconstructed_values[kk][0][jj] = result[0];
          x_reconstructed_values[kk][1][jj] = result[1];
        }
      }

      // reconstruction in y' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      // third index: z' direction
      const auto& num_quad_points = quadrature.size();
      FieldVector<std::vector<FieldVector<RangeType, stencil_size>>, 2> y_reconstructed_values(
          (std::vector<FieldVector<RangeType, stencil_size>>(num_quad_points)));
      result.resize(num_quad_points);
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t ii = 0; ii < 2; ++ii) {
          slope_reconstruction(x_reconstructed_values[kk][ii], result, quadrature);
          for (size_t jj = 0; jj < num_quad_points; ++jj)
            y_reconstructed_values[ii][jj][kk] = result[jj];
        } // ii
      } // kk

      // reconstruction in z' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      // third index: quadrature_points in z' direction
      FieldVector<std::vector<std::vector<RangeType>>, 2> reconstructed_values(
          std::vector<std::vector<RangeType>>(num_quad_points, std::vector<RangeType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < num_quad_points; ++jj)
          slope_reconstruction(y_reconstructed_values[ii][jj], reconstructed_values[ii][jj], quadrature);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < num_quad_points; ++kk) {
            // convert back to non-characteristic variables and to FieldVector instead of EigenVector
            eigenvectors.mv(reconstructed_values[ii][jj][kk], value);
            IntersectionLocalCoordType quadrature_point = {quadrature[jj].position(), quadrature[kk].position()};
            reconstructed_values_map.insert(
                std::make_pair(intersections[2 * dd + ii].geometryInInside().global(quadrature_point), value));
          } // kk
        } // jj
      } // ii
    } // static void reconstruct(...)
  }; // struct helper<3, ...>


  class StencilIterator
  {
  public:
    static constexpr size_t stencil_x = stencil[0];
    static constexpr size_t stencil_y = stencil[1];
    static constexpr size_t stencil_z = stencil[2];

    static void apply(const std::vector<RangeType>& source_values,
                      const BoundaryValueType& boundary_values,
                      FieldVector<FieldVector<FieldVector<RangeType, stencil_z>, stencil_y>, stencil_x>& values,
                      const EntityType& entity,
                      const GridLayerType& grid_layer,
                      const int direction,
                      FieldVector<int, 3> offsets)
    {
      const auto& entity_index = grid_layer.indexSet().index(entity);
      values[stencil_x / 2 + offsets[0]][stencil_y / 2 + offsets[1]][stencil_z / 2 + offsets[2]] =
          source_values[entity_index];
      std::vector<int> boundary_dirs;
      for (const auto& intersection : Dune::intersections(grid_layer, entity)) {
        const auto& intersection_index = intersection.indexInInside();
        if (!end_of_stencil(intersection_index, offsets)) {
          auto new_offsets = offsets;
          if (intersection.boundary() && !intersection.neighbor()) {
            boundary_dirs.push_back(intersection_index);
            const auto& boundary_value =
                boundary_values.local_function(entity)->evaluate(intersection.geometryInInside().global(
                    intersection.geometry().local(intersection.geometry().center())));
            while (!end_of_stencil(intersection_index, new_offsets)) {
              walk(intersection_index, new_offsets);
              values[stencil_x / 2 + new_offsets[0]][stencil_y / 2 + new_offsets[1]][stencil_z / 2 + new_offsets[2]] =
                  boundary_value;
            }
          } else if (direction_allowed(direction, intersection_index)) {
            const auto& outside = intersection.outside();
            walk(intersection_index, new_offsets);
            StencilIterator::apply(
                source_values, boundary_values, values, outside, grid_layer, intersection_index, new_offsets);
          }
        }
      } // intersections

      // TODO: improve multiple boundary handling, currently everything is filled with the boundary value in the first
      // direction
      assert(boundary_dirs.size() <= 3);
      if (boundary_dirs.size() > 1) {
        walk(boundary_dirs[0], offsets);
        const auto& boundary_value =
            values[stencil_x / 2 + offsets[0]][stencil_y / 2 + offsets[1]][stencil_z / 2 + offsets[2]];
        for (auto& values_yz : values)
          for (auto& values_z : values_yz)
            for (auto& value : values_z)
              if (value.size() == 0)
                value = boundary_value;
      }
    } // void apply(...)

  private:
    static void walk(const int dir, FieldVector<int, 3>& offsets)
    {
      dir % 2 ? offsets[dir / 2]++ : offsets[dir / 2]--;
    }

    // Direction is allowed if end of stencil is not reached and direction is not visited by another iterator.
    // Iterators never change direction, they may only spawn new iterators in the directions that have a higher
    // index (i.e. iterators that walk in x direction will spawn iterators going in y and z direction,
    // iterators going in y direction will only spawn iterators in z-direction and z iterators only walk
    // without emitting new iterators).
    static bool direction_allowed(const int dir, const int new_dir)
    {
      return dir == -1 || new_dir == dir || new_dir / 2 > dir / 2;
    }

    static bool end_of_stencil(const int dir, const FieldVector<size_t, 3>& offsets)
    {
      return !(offsets[dir / 2] < stencil[dir / 2] / 2);
    }
  }; // class StencilIterator<...>

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   std::vector<RangeType>& result,
                                   const QuadratureType& quadrature)
  {
    assert(stencil_size == 3);
    std::vector<DomainFieldType> points(quadrature.size());
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      points[ii] = quadrature[ii].position();
    assert(result.size() == points.size());
    const auto& u_left = cell_values[0];
    const auto& u_entity = cell_values[1];
    const auto& u_right = cell_values[2];
    const auto slope_left = u_entity - u_left;
    const auto slope_right = u_right - u_entity;
    auto slope_center = u_right - u_left;
    slope_center *= 0.5;
    const auto slope = internal::ChooseLimiter<slope_limiter>::limit(slope_left, slope_right, slope_center);
    std::fill(result.begin(), result.end(), u_entity);
    for (size_t ii = 0; ii < points.size(); ++ii) {
      auto slope_scaled = slope;
      slope_scaled *= points[ii] - 0.5;
      result[ii] += slope_scaled;
    }
  }

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   FieldVector<RangeType, 2>& result,
                                   const QuadratureType& quadrature)
  {
    std::vector<RangeType> result_vec(2);
    slope_reconstruction(cell_values, result_vec, quadrature);
    std::copy(result_vec.begin(), result_vec.end(), result.begin());
  }

  const std::vector<RangeType> source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const bool is_linear_;
  const QuadratureType quadrature_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values_;
  static thread_local std::unique_ptr<EigenSolverType> eigensolver_;
  static bool is_instantiated_;
}; // class LocalReconstructionFvOperator

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter,
          class EigenSolverType>
constexpr std::array<size_t, 3> LocalReconstructionFvOperator<GridLayerType,
                                                              AnalyticalFluxType,
                                                              BoundaryValueType,
                                                              polOrder,
                                                              slope_limiter,
                                                              EigenSolverType>::stencil;

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter,
          class EigenSolverType>
thread_local std::unique_ptr<EigenSolverType> LocalReconstructionFvOperator<GridLayerType,
                                                                            AnalyticalFluxType,
                                                                            BoundaryValueType,
                                                                            polOrder,
                                                                            slope_limiter,
                                                                            EigenSolverType>::eigensolver_;

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter,
          class EigenSolverType>
bool LocalReconstructionFvOperator<GridLayerType,
                                   AnalyticalFluxType,
                                   BoundaryValueType,
                                   polOrder,
                                   slope_limiter,
                                   EigenSolverType>::is_instantiated_(false);


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH
