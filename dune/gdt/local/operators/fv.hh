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

#ifndef DUNE_GDT_LOCAL_OPERATORS_FV_HH
#define DUNE_GDT_LOCAL_OPERATORS_FV_HH

#include "interfaces.hh"

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/container/eigen.hh>

namespace Dune {
namespace GDT {

enum class SlopeLimiters
{
  minmod,
  mc,
  superbee,
  no_slope
};

// forwards
template <class NumericalFluxType>
class LocalCouplingFvOperator;

template <class NumericalFluxType>
class LocalBoundaryFvOperator;

template <class MatrixImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class LocalReconstructionFvOperator;

template <class GridViewType,
          class AnalyticalFluxType,
          size_t dimDomain,
          size_t dimRange,
          size_t polOrder,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class LocalWENOReconstructionFvOperator;


namespace internal {


// Traits
template <class NumericalFluxType>
struct LocalCouplingFvOperatorTraits
{
  typedef LocalCouplingFvOperator<NumericalFluxType> derived_type;
};

template <class NumericalFluxType>
struct LocalBoundaryFvOperatorTraits
{
  typedef LocalBoundaryFvOperator<NumericalFluxType> derived_type;
};

template <class MatrixImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
struct LocalReconstructionFvOperatorTraits
{
  typedef LocalReconstructionFvOperator<MatrixImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
  typedef MatrixImp MatrixType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::RangeFieldType RangeFieldType;
  static const size_t dimRange = BoundaryValueFunctionType::dimRange;
};

// template<class GridViewType, size_t domainDim, size_t rangeDim, size_t polOrder>
// struct LocalWENOReconstructionFvOperatorTraits
//{
//  static const size_t dimRange = rangeDim;
//  static const size_t dimDomain = domainDim;
//  typedef LocalWENOReconstructionFvOperator<GridViewType, domainDim, rangeDim, polOrder> derived_type;
//  typedef typename GridViewType::ctype RangeFieldType;
//  typedef FieldVector<RangeFieldType, dimRange> RangeType;
//};

struct FieldVectorLess
{
  template <class FieldType, int dimDomain>
  bool operator()(const FieldVector<FieldType, dimDomain>& a, const FieldVector<FieldType, dimDomain>& b) const
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (XT::Common::FloatCmp::lt(a[dd], b[dd]))
        return true;
      else if (XT::Common::FloatCmp::gt(a[dd], b[dd]))
        return false;
    }
    return false;
  }
};

template <SlopeLimiters slope_limiter, class VectorType>
struct ChooseLimiter
{
  static VectorType
  limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& centered_slope);
};

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::minmod, VectorType>
{
  static VectorType
  limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& /*centered_slope*/)
  {
    VectorType ret;
    for (size_t ii = 0; ii < slope_left.size(); ++ii) {
      const auto slope_left_abs = std::abs(slope_left[ii]);
      const auto slope_right_abs = std::abs(slope_right[ii]);
      if (slope_left_abs < slope_right_abs && slope_left[ii] * slope_right[ii] > 0)
        ret[ii] = slope_left[ii];
      else if (Dune::XT::Common::FloatCmp::ge(slope_left_abs, slope_right_abs) && slope_left[ii] * slope_right[ii] > 0)
        ret[ii] = slope_right[ii];
      else
        ret[ii] = 0.0;
    }
    return ret;
  }
};

// template <class FieldType>
// struct ChooseLimiter<SlopeLimiters::minmod, XT::LA::EigenDenseVector<FieldType>>
//{
//  typedef XT::LA::EigenDenseVector<FieldType> VectorType;
//  static VectorType
//  limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& /*centered_slope*/)
//  {
//    VectorType ret(slope_left.size());
//    for (size_t ii = 0; ii < slope_left.size(); ++ii) {
//      const auto slope_left_abs = std::abs(slope_left.get_entry(ii));
//      const auto slope_right_abs = std::abs(slope_right.get_entry(ii));
//      if (XT::Common::FloatCmp::lt(slope_left_abs, slope_right_abs)
//          && slope_left.get_entry(ii) * slope_right.get_entry(ii) > 0)
//        ret.set_entry(ii, slope_left.get_entry(ii));
//      else if (Dune::XT::Common::FloatCmp::gt(slope_left_abs, slope_right_abs)
//               && slope_left.get_entry(ii) * slope_right.get_entry(ii) > 0)
//        ret.set_entry(ii, slope_right.get_entry(ii));
//      else if (Dune::XT::Common::FloatCmp::eq(slope_left_abs, slope_right_abs))
//        ret.set_entry(ii, (slope_right.get_entry(ii) + slope_left.get_entry(ii)) * 0.5);
//      else
//        ret.set_entry(ii, 0.);
//    }
//    return ret;
//  }
//};

template <class FieldType>
struct ChooseLimiter<SlopeLimiters::minmod, XT::LA::EigenDenseVector<FieldType>>
{
  typedef XT::LA::EigenDenseVector<FieldType> VectorType;
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& centered_slope)
  {
    VectorType ret(slope_left.size(), 0);
    for (size_t ii = 0; ii < slope_left.size(); ++ii) {
      // check for equal sign
      if (slope_left.get_entry(ii) * slope_right.get_entry(ii) > 0
          && centered_slope.get_entry(ii) * slope_right.get_entry(ii) > 0) {
        const auto slope_left_abs = std::abs(slope_left.get_entry(ii));
        const auto slope_right_abs = std::abs(slope_right.get_entry(ii));
        const auto slope_centered_abs = std::abs(centered_slope.get_entry(ii));
        if (XT::Common::FloatCmp::lt(slope_left_abs, slope_right_abs)) {
          if (XT::Common::FloatCmp::lt(slope_left_abs, slope_centered_abs))
            ret.set_entry(ii, slope_left.get_entry(ii));
          else
            ret.set_entry(ii, centered_slope.get_entry(ii));
        } else if (Dune::XT::Common::FloatCmp::gt(slope_left_abs, slope_right_abs)) {
          if (XT::Common::FloatCmp::lt(slope_right_abs, slope_centered_abs))
            ret.set_entry(ii, slope_right.get_entry(ii));
          else
            ret.set_entry(ii, centered_slope.get_entry(ii));
        } else {
        }
      }
    }
    return ret;
  }
};


template <class VectorType>
struct ChooseLimiter<SlopeLimiters::superbee, VectorType>
{
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& centered_slope)
  {
    typedef ChooseLimiter<SlopeLimiters::minmod, VectorType> MinmodType;
    return maxmod(MinmodType::limit(slope_left, slope_right * 2.0, centered_slope),
                  MinmodType::limit(slope_left * 2.0, slope_right, centered_slope));
  }

  static VectorType maxmod(const VectorType& slope_left, const VectorType& slope_right)
  {
    VectorType ret;
    for (size_t ii = 0; ii < slope_left.size(); ++ii) {
      const auto slope_left_abs = std::abs(slope_left[ii]);
      const auto slope_right_abs = std::abs(slope_right[ii]);
      if (slope_left_abs > slope_right_abs && slope_left[ii] * slope_right[ii] > 0)
        ret[ii] = slope_left[ii];
      else if (Dune::XT::Common::FloatCmp::le(slope_left_abs, slope_right_abs) && slope_left[ii] * slope_right[ii] > 0)
        ret[ii] = slope_right[ii];
      else
        ret[ii] = 0.0;
    }
    return ret;
  }
};

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::mc, VectorType>
{
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& centered_slope)
  {
    typedef ChooseLimiter<SlopeLimiters::minmod, VectorType> MinmodType;
    return MinmodType::limit(
        MinmodType::limit(slope_left * 2.0, slope_right * 2.0, centered_slope), centered_slope, centered_slope);
  }
};

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::no_slope, VectorType>
{
  static VectorType
  limit(const VectorType& /*slope_left*/, const VectorType& /*slope_right*/, const VectorType& /*centered_slope*/)
  {
    return VectorType(0.);
  }
};

template <class FieldType>
struct ChooseLimiter<SlopeLimiters::no_slope, XT::LA::EigenDenseVector<FieldType>>
{
  typedef XT::LA::EigenDenseVector<FieldType> VectorType;
  static VectorType
  limit(const VectorType& slope_left, const VectorType& /*slope_right*/, const VectorType& /*centered_slope*/)
  {
    return VectorType(slope_left.size(), 0.);
  }
};


} // namespace internal


template <class NumericalFluxType>
class LocalCouplingFvOperator
    : public LocalCouplingOperatorInterface<internal::LocalCouplingFvOperatorTraits<NumericalFluxType>>
{
public:
  typedef typename NumericalFluxType::RangeFieldType RangeFieldType;
  static const size_t dimDomain = NumericalFluxType::dimDomain;
  static const size_t dimRange = NumericalFluxType::dimRange;
  typedef Dune::QuadratureRule<RangeFieldType, dimDomain - 1> QuadratureRuleType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> ResultType;


  template <class... Args>
  explicit LocalCouplingFvOperator(const QuadratureRuleType& quadrature_rule, Args&&... args)
    : quadrature_rule_(quadrature_rule)
    , numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor) const
  {
    const auto entity = intersection.inside();
    const auto neighbor = intersection.outside();
    const auto local_source_entity = source.local_function(entity);
    const auto local_source_neighbor = source.local_function(neighbor);
    const auto local_functions_tuple_entity = numerical_flux_.local_functions(entity);
    const auto local_functions_tuple_neighbor = numerical_flux_.local_functions(neighbor);
    ResultType result(0);
    for (const auto& quad_point : quadrature_rule_) {
      result += ResultType(numerical_flux_.evaluate(local_functions_tuple_entity,
                                                    local_functions_tuple_neighbor,
                                                    *local_source_entity,
                                                    *local_source_neighbor,
                                                    intersection,
                                                    quad_point.position()))
                * quad_point.weight();
    }
    local_range_entity.vector().add(result * (1.0 / entity.geometry().volume()));
    local_range_neighbor.vector().add(result * (-1.0 / neighbor.geometry().volume()));
  }

private:
  const QuadratureRuleType& quadrature_rule_;
  const NumericalFluxType numerical_flux_;
};


template <class NumericalFluxType>
class LocalBoundaryFvOperator
    : public LocalBoundaryOperatorInterface<internal::LocalBoundaryFvOperatorTraits<NumericalFluxType>>
{
public:
  typedef typename NumericalFluxType::RangeFieldType RangeFieldType;
  static const size_t dimDomain = NumericalFluxType::dimDomain;
  static const size_t dimRange = NumericalFluxType::dimRange;
  typedef Dune::QuadratureRule<RangeFieldType, dimDomain - 1> QuadratureRuleType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> ResultType;

  template <class... Args>
  explicit LocalBoundaryFvOperator(const QuadratureRuleType& quadrature_rule, Args&&... args)
    : quadrature_rule_(quadrature_rule)
    , numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity) const
  {
    const auto entity = intersection.inside();
    const auto local_source_entity = source.local_function(entity);
    const auto local_functions_tuple = numerical_flux_.local_functions(entity);
    typedef
        typename Dune::XT::Common::FieldVector<typename LocalDiscreteFunction<SpaceType, VectorType>::DomainFieldType,
                                               LocalDiscreteFunction<SpaceType, VectorType>::dimRange>
            ResultType;
    ResultType result(0);
    for (const auto& quad_point : quadrature_rule_) {
      result += ResultType(numerical_flux_.evaluate(
                    local_functions_tuple, *local_source_entity, intersection, quad_point.position()))
                * quad_point.weight();
    }
    local_range_entity.vector().add(result / entity.geometry().volume());
  }

private:
  const QuadratureRuleType& quadrature_rule_;
  const NumericalFluxType numerical_flux_;
};

enum class BasisFunctionType
{
  legendre,
  hat_functions,
  first_order_dg
};

template <class RangeType, BasisFunctionType basis_function_type>
struct rescale_factor
{
  static typename RangeType::field_type get(const RangeType& /*u*/, const RangeType& /*u_bar*/)
  {
    DUNE_THROW(NotImplemented, "Not implemented for this basisfunction type");
    return 0;
  }
};

template <class RangeType>
class rescale_factor<RangeType, BasisFunctionType::legendre>
{
  static typename RangeType::field_type get(const RangeType& u, const RangeType& u_bar)
  {
    return 2 * std::max(u[0], u_bar[0]);
  }
};

template <class RangeType>
struct rescale_factor<RangeType, BasisFunctionType::hat_functions>
{
  static typename RangeType::field_type get(const RangeType& u, const RangeType& u_bar)
  {
    typedef typename RangeType::field_type FieldType;
    return 2 * std::max(std::accumulate(u.begin(), u.end(), FieldType(0)),
                        std::accumulate(u_bar.begin(), u_bar.end(), FieldType(0)));
  }
};

template <class RangeType>
struct rescale_factor<RangeType, BasisFunctionType::first_order_dg>
{
  static typename RangeType::field_type get(const RangeType& u, const RangeType& u_bar)
  {
    typename RangeType::FieldType u_sum;
    auto u_bar_sum = u_sum;
    for (size_t ii = 0; ii < u.size(); ii += 4) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
  }
};


template <class GridViewType,
          size_t dimDomain,
          size_t dimRange,
          BasisFunctionType basis_function_type = BasisFunctionType::legendre>
class LocalRealizabilityLimiter : public XT::Grid::Functor::Codim0<GridViewType>
{
  typedef typename GridViewType::ctype FieldType;
  typedef typename XT::LA::EigenDenseVector<FieldType> EigenVectorType;
  typedef FieldVector<FieldType, dimRange> RangeType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit LocalRealizabilityLimiter(
      const GridViewType& grid_view,
      const std::vector<FieldVector<FieldType, dimRange + 1>>& plane_coefficients,
      const std::vector<std::vector<std::vector<EigenVectorType>>>& cell_averages,
      const std::vector<FieldVector<size_t, dimDomain>>& entity_indices,
      std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                           RangeType,
                           internal::FieldVectorLess>>& reconstructed_values)
    : grid_view_(grid_view)
    , plane_coefficients_(plane_coefficients)
    , cell_averages_(cell_averages)
    , entity_indices_(entity_indices)
    , reconstructed_values_(reconstructed_values)
  {
  }

  void apply_local(const EntityType& entity)
  {
    const auto& entity_index = grid_view_.indexSet().index(entity);
    auto& local_reconstructed_values = reconstructed_values_[entity_index];

    // get cell average
    const size_t i_x = entity_indices_[entity_index][0];
    const size_t i_y = entity_indices_[entity_index][1];
    const size_t i_z = entity_indices_[entity_index][2];
    const auto& u_bar_eigen = cell_averages_[i_x][i_y][i_z];
    assert(u_bar_eigen.size() == dimRange);
    RangeType u_bar;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u_bar[ii] = u_bar_eigen.get_entry(ii);

    const FieldType epsilon = 1e-3;

    // vector to store thetas for each local reconstructed value
    std::vector<FieldType> thetas(local_reconstructed_values.size(), -epsilon);

    FieldVector<FieldType, dimRange> a;
    FieldType b;
    for (const auto& coeffs : plane_coefficients_) {
      std::copy_n(coeffs.begin(), dimRange, a.begin());
      b = coeffs[dimRange];
      size_t ll = -1;
      for (const auto& pair : local_reconstructed_values) {
        ++ll;
        auto u_l = pair.second;
        // if u_l equals u_bar, no limiting is necessary
        if (XT::Common::FloatCmp::eq(u_l, u_bar))
          continue;
        // rescale u_l, u_bar
        const auto factor = rescale_factor<RangeType, basis_function_type>::get(u_l, u_bar);
        u_l /= factor;
        auto u_bar_l = u_bar;
        u_bar_l /= factor;

        FieldType theta_li = (b - a * u_l) / (a * (u_bar_l - u_l));
        if (XT::Common::FloatCmp::ge(theta_li, -epsilon) && XT::Common::FloatCmp::le(theta_li, 1.))
          thetas[ll] = std::max(thetas[ll], theta_li);
      }
    }
    for (auto& theta : thetas)
      theta = std::min(epsilon + theta, 1.);

    auto theta_entity = *std::max(thetas.begin(), thetas.end());

    for (auto& pair : local_reconstructed_values) {
      auto& u = pair.second;
      auto u_scaled = u;
      u_scaled *= (1 - theta_entity);
      auto u_bar_scaled = u_bar;
      u_bar_scaled *= theta_entity;
      u = u_scaled + u_bar_scaled;
    }
  } // void apply_local(...)

private:
  const GridViewType& grid_view_;
  const std::vector<FieldVector<FieldType, dimRange + 1>>& plane_coefficients_;
  const std::vector<std::vector<std::vector<EigenVectorType>>>& cell_averages_;
  const std::vector<FieldVector<size_t, 3>>& entity_indices_;
  std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                       RangeType,
                       internal::FieldVectorLess>>& reconstructed_values_;
}; // class LocalRealizabilityLimiter

template <class GridViewType,
          class AnalyticalFluxType,
          size_t dimDomain,
          size_t dimRange,
          size_t polOrder,
          SlopeLimiters slope_limiter>
class LocalWENOReconstructionFvOperator : public XT::Grid::Functor::Codim0<GridViewType>
{
  typedef typename GridViewType::ctype FieldType;
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static const size_t stencil_size = 2 * polOrder + 1;
  typedef typename Dune::XT::Common::FieldVector<FieldType, dimRange> XTFieldVectorType;
  typedef typename XT::LA::EigenDenseVector<FieldType> EigenVectorType;
  typedef typename XT::LA::EigenDenseMatrix<FieldType> EigenMatrixType;
  typedef FieldVector<FieldType, dimRange> RangeType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit LocalWENOReconstructionFvOperator(
      const GridViewType& grid_view,
      const AnalyticalFluxType& analytical_flux,
      const XT::Common::Parameter& param,
      const std::vector<std::vector<std::vector<EigenVectorType>>>& cell_averages,
      const std::vector<FieldVector<size_t, dimDomain>>& entity_indices,
      const FieldVector<Dune::QuadratureRule<FieldType, 1>, dimDomain>& quadrature_rules,
      std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                           RangeType,
                           internal::FieldVectorLess>>& reconstructed_values)
    : grid_view_(grid_view)
    , analytical_flux_(analytical_flux)
    , t_(param.get("t")[0])
    , cell_averages_(cell_averages)
    , entity_indices_(entity_indices)
    , quadrature_rules_(quadrature_rules)
    , reconstructed_values_(reconstructed_values)
  {
  }

  template <int stencil_size>
  void slope_reconstruction(const FieldVector<EigenVectorType, stencil_size>& cell_values,
                            std::vector<EigenVectorType>& result,
                            Dune::QuadratureRule<FieldType, 1> quadrature_rule)
  {
    static_assert(stencil_size == 3, "");
    std::vector<FieldType> points(quadrature_rule.size());
    for (size_t ii = 0; ii < quadrature_rule.size(); ++ii)
      points[ii] = quadrature_rule[ii].position();
    assert(result.size() == points.size());
    const auto& u_left = cell_values[0];
    const auto& u_entity = cell_values[1];
    const auto& u_right = cell_values[2];
    const auto slope_left = u_entity - u_left;
    const auto slope_right = u_right - u_entity;
    const auto slope_centered = (u_right - u_left) * 0.5;
    const auto slope =
        internal::ChooseLimiter<slope_limiter, EigenVectorType>::limit(slope_left, slope_right, slope_centered);
    for (size_t ii = 0; ii < points.size(); ++ii)
      result[ii] = u_entity + slope * (points[ii] - 0.5);
  }

  void apply_local(const EntityType& entity)
  {
    const auto& entity_center = entity.geometry().center();
    const auto& entity_index = grid_view_.indexSet().index(entity);
    auto& local_reconstructed_values = reconstructed_values_[entity_index];
    // get cell averages on stencil
    const size_t i_x = entity_indices_[entity_index][0];
    const size_t i_y = entity_indices_[entity_index][1];
    const size_t i_z = entity_indices_[entity_index][2];
    static const size_t offset_0 = polOrder;
    static const size_t offset_1 = polOrder;
    static const size_t offset_2 = polOrder;
    FieldVector<FieldVector<FieldVector<EigenVectorType, stencil_size>, stencil_size>, stencil_size> values;
    const FieldVector<size_t, dimDomain> grid_sizes{
        cell_averages_.size() - 2, cell_averages_[0].size() - 2, cell_averages_[0][0].size() - 2};
    for (size_t ii = 0; ii < stencil_size; ++ii) {
      const size_t index_x = i_x + ii >= offset_0
                                 ? (i_x - offset_0 + ii < grid_sizes[0] ? i_x - offset_0 + ii : grid_sizes[0] + 1)
                                 : grid_sizes[0];
      auto& values_x = values[ii];
      const auto& averages_x = cell_averages_[index_x];
      for (size_t jj = 0; jj < stencil_size; ++jj) {
        const size_t index_y = i_y + jj >= offset_1
                                   ? (i_y - offset_1 + jj < grid_sizes[1] ? i_y - offset_1 + jj : grid_sizes[1] + 1)
                                   : grid_sizes[1];
        auto& values_xy = values_x[jj];
        const auto& averages_xy = averages_x[index_y];
        for (size_t kk = 0; kk < stencil_size; ++kk) {
          const size_t index_z = i_z + kk >= offset_2
                                     ? (i_z - offset_2 + kk < grid_sizes[2] ? i_z - offset_2 + kk : grid_sizes[2] + 1)
                                     : grid_sizes[2];
          values_xy[kk] = averages_xy[index_z];
        } // kk (z-dimension)
      } // jj (y-dimension)
    } // ii (x-dimension)

    // get jacobians
    const auto& u_entity_eigen = cell_averages_[i_x][i_y][i_z];
    FieldVector<FieldType, dimRange> u_entity;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u_entity[ii] = u_entity_eigen.get_entry(ii);
    const FieldVector<FieldMatrix<FieldType, dimRange, dimRange>, dimDomain> jacobians =
        analytical_flux_.jacobian(u_entity, entity, entity.geometry().local(entity_center), t_);

    // get intersections
    FieldVector<typename GridViewType::Intersection, 2 * dimDomain> faces;
    for (const auto& intersection : Dune::intersections(grid_view_, entity)) {
      const auto& n = intersection.unitOuterNormal(intersection.geometry().local(intersection.geometry().center()));
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        if (XT::Common::FloatCmp::eq(n[dd], -1.))
          faces[2 * dd] = intersection;
        else if (XT::Common::FloatCmp::eq(n[dd], 1.))
          faces[2 * dd + 1] = intersection;
      }
    }

    // do WENO reconstruction
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // get jacobian and convert to EigenMatrixType;
      // TODO: proper conversion without strings
      const auto& jacobian = jacobians[dd];
      static const size_t precision = 20;
      const auto jacobian_eigen =
          Dune::XT::Common::from_string<EigenMatrixType>(Dune::XT::Common::to_string(jacobian, precision));

      // get eigenvectors
      ::Eigen::SelfAdjointEigenSolver<typename EigenMatrixType::BackendType> eigen_solver(jacobian_eigen.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvectors_backend = eigen_solver.eigenvectors();
      const auto eigenvectors = EigenMatrixType(eigenvectors_backend.real());
      const auto eigenvectors_inverse = EigenMatrixType(eigenvectors_backend.inverse().real());

      // set x-direction to x[dd] and set y- and z-direction
      const Dune::QuadratureRule<FieldType, 1>& quadrature_y = quadrature_rules_[(dd + 1) % 3];
      const Dune::QuadratureRule<FieldType, 1>& quadrature_z = quadrature_rules_[(dd + 2) % 3];

      // convert to characteristic variables
      FieldVector<FieldVector<FieldVector<EigenVectorType, stencil_size>, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        for (size_t jj = 0; jj < stencil_size; ++jj)
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            if (dd == 0)
              char_values[ii][jj][kk] = eigenvectors_inverse * values[ii][jj][kk];
            else if (dd == 1)
              char_values[jj][kk][ii] = eigenvectors_inverse * values[ii][jj][kk];
            else if (dd == 2)
              char_values[kk][ii][jj] = eigenvectors_inverse * values[ii][jj][kk];
            else
              DUNE_THROW(NotImplemented, "Not implemented in more than 3 dimensions!");
          }

      // reconstruction in x direction
      // first index: z direction
      // second index: dimension 2 for left/right interface
      // third index: y direction
      FieldVector<FieldVector<FieldVector<EigenVectorType, stencil_size>, 2>, stencil_size> first_reconstructed_values;
      std::vector<EigenVectorType> result(2);
      // quadrature rule containing left and right interface points
      Dune::QuadratureRule<FieldType, 1> points;
      points.push_back(Dune::QuadraturePoint<FieldType, 1>(0., 0.5));
      points.push_back(Dune::QuadraturePoint<FieldType, 1>(1., 0.5));
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          //         WENOreconstruction(char_values[ii][jj], result);
          slope_reconstruction(char_values[ii][jj], result, points);
          first_reconstructed_values[ii][0][jj] = result[0];
          first_reconstructed_values[ii][1][jj] = result[1];
        }
      }

      const size_t num_quad_points_y = quadrature_y.size();
      // reconstruction in y direction
      // first index: left/right interface
      // second index: quadrature_points in y direction
      // third index: z direction
      FieldVector<std::vector<FieldVector<EigenVectorType, stencil_size>>, 2> second_reconstructed_values(
          (std::vector<FieldVector<EigenVectorType, stencil_size>>(num_quad_points_y)));
      result.resize(num_quad_points_y);
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t kk = 0; kk < 2; ++kk) {
          //          WENOreconstruction(first_reordered[ii][kk], result, 1d_quadrature_points_y);
          slope_reconstruction(first_reconstructed_values[ii][kk], result, quadrature_y);
          for (size_t ll = 0; ll < num_quad_points_y; ++ll)
            second_reconstructed_values[kk][ll][ii] = result[ll];
        } // kk
      } // ii

      const size_t num_quad_points_z = quadrature_z.size();
      // reconstruction in z direction
      // first index: left/right interface
      // second index: quadrature_points in y direction
      // third index: quadrature_points in z direction
      FieldVector<std::vector<std::vector<EigenVectorType>>, 2> reconstructed_values(
          std::vector<std::vector<EigenVectorType>>(num_quad_points_y,
                                                    std::vector<EigenVectorType>(num_quad_points_z)));
      for (size_t kk = 0; kk < 2; ++kk) {
        for (size_t ll = 0; ll < num_quad_points_y; ++ll) {
          //          WENOreconstruction(second_reconstructed_values[kk][ll], reconstructed_values[kk][ll],
          //          1d_quadrature_points_z);
          slope_reconstruction(second_reconstructed_values[kk][ll], reconstructed_values[kk][ll], quadrature_z);
        }
      }

      // convert coordinates on face to local entity coordinates and store
      typedef typename GridViewType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;
      for (size_t kk = 0; kk < 2; ++kk) {
        for (size_t ll = 0; ll < num_quad_points_y; ++ll) {
          for (size_t mm = 0; mm < num_quad_points_z; ++mm) {
            // convert back to non-characteristic variables and to FieldVector instead of EigenVector
            const auto value = XT::Common::from_string<RangeType>(
                XT::Common::to_string(eigenvectors * reconstructed_values[kk][ll][mm], precision));
            IntersectionLocalCoordType quadrature_point{quadrature_y[ll].position()[0], quadrature_z[mm].position()[0]};
            local_reconstructed_values.insert(
                std::make_pair(faces[2 * dd + kk].geometryInInside().global(quadrature_point), value));
          } // mm
        } // ll
      } // kk
    } // dd
  } // void apply_local(...)

private:
  const GridViewType& grid_view_;
  const AnalyticalFluxType& analytical_flux_;
  const double t_;
  const std::vector<std::vector<std::vector<EigenVectorType>>>& cell_averages_;
  const std::vector<FieldVector<size_t, 3>>& entity_indices_;
  const FieldVector<Dune::QuadratureRule<FieldType, 1>, dimDomain>& quadrature_rules_;
  std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                       RangeType,
                       internal::FieldVectorLess>>& reconstructed_values_;
}; // class LocalWENOReconstructionFvOperator

template <class GridViewType, class AnalyticalFluxType, size_t dimRange, size_t polOrder, SlopeLimiters slope_limiter>
class LocalWENOReconstructionFvOperator<GridViewType, AnalyticalFluxType, 1, dimRange, polOrder, slope_limiter>
    : public XT::Grid::Functor::Codim0<GridViewType>
{
  static const size_t dimDomain = 1;
  typedef typename GridViewType::ctype FieldType;
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static const size_t stencil_size = 2 * polOrder + 1;
  typedef typename XT::LA::EigenDenseVector<FieldType> EigenVectorType;
  typedef typename XT::LA::EigenDenseMatrix<FieldType> EigenMatrixType;
  typedef FieldVector<FieldType, dimRange> RangeType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit LocalWENOReconstructionFvOperator(
      const GridViewType& grid_view,
      const AnalyticalFluxType& analytical_flux,
      const XT::Common::Parameter& param,
      const std::vector<EigenVectorType>& cell_averages,
      const std::vector<FieldVector<size_t, 1>>& entity_indices,
      const FieldVector<Dune::QuadratureRule<FieldType, 1>, dimDomain>& /*quadrature_rules*/,
      std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                           RangeType,
                           internal::FieldVectorLess>>& reconstructed_values)
    : grid_view_(grid_view)
    , analytical_flux_(analytical_flux)
    , t_(param.get("t")[0])
    , cell_averages_(cell_averages)
    , entity_indices_(entity_indices)
    , reconstructed_values_(reconstructed_values)
  {
  }

  template <int stencil_size>
  void slope_reconstruction(const FieldVector<EigenVectorType, stencil_size>& cell_values,
                            std::vector<EigenVectorType>& result)
  {
    static_assert(stencil_size == 3, "");
    FieldVector<FieldType, 2> points{0., 1.};
    assert(result.size() == points.size());

    const auto& u_left = cell_values[0];
    const auto& u_entity = cell_values[1];
    const auto& u_right = cell_values[2];
    const auto slope_left = u_entity - u_left;
    const auto slope_right = u_right - u_entity;
    const auto slope_centered = (u_right - u_left) * 0.5;
    const auto slope =
        internal::ChooseLimiter<slope_limiter, EigenVectorType>::limit(slope_left, slope_right, slope_centered);
    for (size_t ii = 0; ii < points.size(); ++ii)
      result[ii] = u_entity + slope * (points[ii] - 0.5);
  }

  void apply_local(const EntityType& entity)
  {
    const auto& entity_index = grid_view_.indexSet().index(entity);
    auto& local_reconstructed_values = reconstructed_values_[entity_index];
    // get cell averages on stencil
    const size_t i_x = entity_indices_[entity_index];
    FieldVector<EigenVectorType, stencil_size> values;
    static const size_t offset = polOrder;
    const size_t grid_size = cell_averages_.size() - 2;
    for (size_t ii = 0; ii < stencil_size; ++ii) {
      const size_t index_x =
          i_x + ii >= offset ? (i_x - offset + ii < grid_size ? i_x - offset + ii : grid_size + 1) : grid_size;
      values[ii] = cell_averages_[index_x];
    } // ii (x-dimension)

    // get jacobians
    const auto& u_entity_eigen = cell_averages_[i_x];
    RangeType u_entity;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u_entity[ii] = u_entity_eigen.get_entry(ii);
    const FieldMatrix<FieldType, dimRange, dimRange> jacobian =
        analytical_flux_.jacobian(u_entity, entity, entity.geometry().local(entity.geometry().center()), t_);

    // get intersections
    std::vector<typename GridViewType::Intersection> faces(2);
    for (const auto& intersection : Dune::intersections(grid_view_, entity)) {
      const auto& n = intersection.unitOuterNormal(intersection.geometry().local(intersection.geometry().center()));
      if (XT::Common::FloatCmp::eq(n[0], -1.))
        faces[0] = intersection;
      else if (XT::Common::FloatCmp::eq(n[0], 1.))
        faces[1] = intersection;
    }

    // do WENO reconstruction
    // get jacobian and convert to EigenMatrixType;
    // TODO: proper conversion without strings
    static const size_t precision = 15;
    const auto jacobian_eigen =
        Dune::XT::Common::from_string<EigenMatrixType>(Dune::XT::Common::to_string(jacobian, precision));

    // get eigenvectors
    ::Eigen::EigenSolver<typename EigenMatrixType::BackendType> eigen_solver(jacobian_eigen.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvectors_backend = eigen_solver.eigenvectors();
#ifndef NDEBUG
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        assert(eigenvectors_backend(ii, jj).imag() < 1e-15);
#endif
    const auto eigenvectors = EigenMatrixType(eigenvectors_backend.real());
    const auto eigenvectors_inverse = EigenMatrixType(eigenvectors_backend.inverse().real());

    // convert to characteristic variables
    auto char_values = values;
    for (auto& value : char_values)
      value = eigenvectors_inverse * value;

    // reconstruction
    std::vector<EigenVectorType> reconstructed_values(2);
    //          WENOreconstruction(char_values[ii][jj], result);
    slope_reconstruction(char_values, reconstructed_values);

    // convert coordinates on face to local entity coordinates and store
    typedef typename GridViewType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;
    for (size_t kk = 0; kk < 2; ++kk) {
      // convert back to non-characteristic variables and to FieldVector instead of EigenVector
      const auto value =
          XT::Common::from_string<RangeType>(XT::Common::to_string(eigenvectors * reconstructed_values[kk], precision));
      local_reconstructed_values.insert(
          std::make_pair(faces[kk].geometryInInside().global(IntersectionLocalCoordType{}), value));
    } // kk
  } // void apply_local(...)

private:
  const GridViewType& grid_view_;
  const AnalyticalFluxType& analytical_flux_;
  const double t_;
  const std::vector<EigenVectorType>& cell_averages_;
  const std::vector<FieldVector<size_t, 1>>& entity_indices_;
  std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                       RangeType,
                       internal::FieldVectorLess>>& reconstructed_values_;
}; // class LocalWENOReconstructionFvOperator

template <class MatrixImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class LocalReconstructionFvOperator
    : public LocalOperatorInterface<internal::LocalReconstructionFvOperatorTraits<MatrixImp,
                                                                                  BoundaryValueFunctionImp,
                                                                                  slope_limiter>>
{
  typedef internal::LocalReconstructionFvOperatorTraits<MatrixImp, BoundaryValueFunctionImp, slope_limiter> Traits;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const size_t dimRange = Traits::dimRange;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> XTFieldVectorType;

public:
  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;

  explicit LocalReconstructionFvOperator(const MatrixType& eigenvectors,
                                         const MatrixType& eigenvectors_inverse,
                                         const std::shared_ptr<BoundaryValueFunctionType>& boundary_values)
    : eigenvectors_(eigenvectors)
    , eigenvectors_inverse_(eigenvectors_inverse)
    , boundary_values_(boundary_values)
  {
  }

  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    const auto& entity = local_range.entity();
    const auto& grid_layer = local_range.space().grid_layer();
    // get reconstruction vector and mapper
    auto& reconstruction_vector = local_range.vector();
    const auto& reconstruction_mapper = local_range.space().mapper();
    const auto entity_center = entity.geometry().center();
    // walk over intersections to get values of discrete function on left and right neighbor entity
    typename SourceType::RangeType u_left, u_right;
    typename SourceType::RangeType u_entity =
        source.local_discrete_function(entity)->evaluate(entity.geometry().local(entity_center));
    const auto i_it_end = grid_layer.iend(entity);
    for (auto i_it = grid_layer.ibegin(entity); i_it != i_it_end; ++i_it) {
      const auto& intersection = *i_it;
      if (intersection.neighbor()) {
        const auto neighbor = intersection.outside();
        const auto neighbor_center = neighbor.geometry().center();
        const bool boundary = intersection.boundary();
        if ((neighbor_center[0] < entity_center[0] && !boundary) || (neighbor_center[0] > entity_center[0] && boundary))
          u_left = source.local_discrete_function(neighbor)->evaluate(neighbor.geometry().local(neighbor_center));
        else
          u_right = source.local_discrete_function(neighbor)->evaluate(neighbor.geometry().local(neighbor_center));
      } else {
        if (intersection.geometry().center()[0] < entity_center[0])
          u_left = boundary_values_->local_function(entity)->evaluate(intersection.geometryInInside().center());
        else
          u_right = boundary_values_->local_function(entity)->evaluate(intersection.geometryInInside().center());
      }
    }

    // diagonalize the system of equations from u_t + A*u_x = 0 to w_t + D*w_x = 0 where D = R^(-1)*A*R, w = R^(-1)*u
    // and R matrix of eigenvectors of A
    const XTFieldVectorType w_left(eigenvectors_inverse_ * u_left);
    const XTFieldVectorType w_right(eigenvectors_inverse_ * u_right);
    const XTFieldVectorType w_entity(eigenvectors_inverse_ * u_entity);

    const XTFieldVectorType w_slope_left = w_entity - w_left;
    const XTFieldVectorType w_slope_right = w_right - w_entity;
    const XTFieldVectorType w_centered_slope = w_right * RangeFieldType(0.5) - w_left * RangeFieldType(0.5);
    const XTFieldVectorType w_slope =
        internal::ChooseLimiter<slope_limiter, XTFieldVectorType>::limit(w_slope_left, w_slope_right, w_centered_slope);
    const XTFieldVectorType half_w_slope = w_slope * RangeFieldType(0.5);
    const XTFieldVectorType w_reconstructed_left = w_entity - half_w_slope;
    const XTFieldVectorType w_reconstructed_right = w_entity + half_w_slope;

    // convert back to u variable
    const XTFieldVectorType u_reconstructed_left(eigenvectors_ * w_reconstructed_left);
    const XTFieldVectorType u_reconstructed_right(eigenvectors_ * w_reconstructed_right);

    for (size_t factor_index = 0; factor_index < dimRange; ++factor_index) {
      // set values on dofs, dof with local index 0 for each factor space corresponds to basis function 1 - x, local
      // index 1 to x
      reconstruction_vector.set(reconstruction_mapper.mapToLocal(factor_index, entity, 0),
                                u_reconstructed_left[factor_index]);
      reconstruction_vector.set(reconstruction_mapper.mapToLocal(factor_index, entity, 1),
                                u_reconstructed_right[factor_index]);
    }
  } // void apply(...)

private:
  const MatrixType& eigenvectors_;
  const MatrixType& eigenvectors_inverse_;
  const std::shared_ptr<BoundaryValueFunctionImp>& boundary_values_;
}; // class LocalReconstructionFvOperator

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_FV_HH
