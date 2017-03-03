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

template <class RHSEvaluationImp>
class LocalRhsFvOperator;

template <class MatrixImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class LocalReconstructionFvOperator;


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

template <class RHSEvaluationImp>
struct LocalRhsFvOperatorTraits
{
  typedef LocalRhsFvOperator<RHSEvaluationImp> derived_type;
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
    return VectorType(0);
  }
};


} // namespace internal


template <class NumericalFluxType>
class LocalCouplingFvOperator
    : public LocalCouplingOperatorInterface<internal::LocalCouplingFvOperatorTraits<NumericalFluxType>>
{
public:
  template <class... Args>
  explicit LocalCouplingFvOperator(Args&&... args)
    : numerical_flux_(std::forward<Args>(args)...)
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
    const auto geometry_intersection = intersection.geometry();
    const auto local_functions_tuple_entity = numerical_flux_.local_functions(entity);
    const auto local_functions_tuple_neighbor = numerical_flux_.local_functions(neighbor);
    const Dune::XT::Common::FieldVector<typename LocalDiscreteFunction<SpaceType, VectorType>::DomainFieldType,
                                        LocalDiscreteFunction<SpaceType, VectorType>::dimRange>
        result = numerical_flux_.evaluate(local_functions_tuple_entity,
                                          local_functions_tuple_neighbor,
                                          *local_source_entity,
                                          *local_source_neighbor,
                                          intersection,
                                          geometry_intersection.local(geometry_intersection.center()));
    local_range_entity.vector().add(result * (1.0 / entity.geometry().volume()));
    local_range_neighbor.vector().add(result * (-1.0 / neighbor.geometry().volume()));
  }

private:
  const NumericalFluxType numerical_flux_;
};


template <class NumericalFluxType>
class LocalBoundaryFvOperator
    : public LocalBoundaryOperatorInterface<internal::LocalBoundaryFvOperatorTraits<NumericalFluxType>>
{
public:
  template <class... Args>
  explicit LocalBoundaryFvOperator(Args&&... args)
    : numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity) const
  {
    const auto entity = intersection.inside();
    const auto local_source_entity = source.local_function(entity);
    const auto geometry_intersection = intersection.geometry();
    const auto local_functions_tuple = numerical_flux_.local_functions(entity);
    auto result = numerical_flux_.evaluate(local_functions_tuple,
                                           *local_source_entity,
                                           intersection,
                                           geometry_intersection.local(geometry_intersection.center()));
    result /= entity.geometry().volume();
    local_range_entity.vector().add(result);
  }

private:
  const NumericalFluxType numerical_flux_;
};

/** TODO: add support for time-dependent RHS
 *  TODO: implement as integral operator??
 * */
template <class RHSEvaluationImp>
class LocalRhsFvOperator : public LocalOperatorInterface<internal::LocalRhsFvOperatorTraits<RHSEvaluationImp>>
{
public:
  explicit LocalRhsFvOperator(const RHSEvaluationImp& rhs_evaluation)
    : rhs_evaluation_(rhs_evaluation)
  {
  }

  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    const auto& entity = local_range.entity();
    const auto local_source_entity = source.local_function(entity);
    const auto x_local = entity.geometry().local(entity.geometry().center());
    const auto u = local_source_entity->evaluate(x_local);
    const auto result = rhs_evaluation_.evaluate(u, entity, x_local);
    local_range.vector().add(result);
  }

private:
  const RHSEvaluationImp& rhs_evaluation_;
};


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
                                         const BoundaryValueFunctionType& boundary_values)
    : eigenvectors_(eigenvectors)
    , eigenvectors_inverse_(eigenvectors_inverse)
    , boundary_values_(boundary_values)
  {
  }

  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    const auto& entity = local_range.entity();
    const auto& grid_view = local_range.space().grid_view();
    // get reconstruction vector and mapper
    auto& reconstruction_vector = local_range.vector();
    const auto& reconstruction_mapper = local_range.space().mapper();
    const auto entity_center = entity.geometry().center();
    // walk over intersections to get values of discrete function on left and right neighbor entity
    typename SourceType::RangeType u_left, u_right;
    typename SourceType::RangeType u_entity =
        source.local_discrete_function(entity)->evaluate(entity.geometry().local(entity_center));
    const auto i_it_end = grid_view.iend(entity);
    for (auto i_it = grid_view.ibegin(entity); i_it != i_it_end; ++i_it) {
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
          u_left = boundary_values_.local_function(entity)->evaluate(intersection.geometryInInside().center());
        else
          u_right = boundary_values_.local_function(entity)->evaluate(intersection.geometryInInside().center());
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
  const BoundaryValueFunctionImp& boundary_values_;
};

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_FV_HH
