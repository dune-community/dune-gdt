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

#include <dune/common/fvector.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/container/eigen.hh>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
//-#include <libqhullcpp/RboxPoints.h>
//-#include <libqhullcpp/QhullError.h>
//-#include <libqhullcpp/QhullQh.h>
//-#include <libqhullcpp/QhullFacet.h>
//-#include <libqhullcpp/QhullLinkedList.h>
//-#include <libqhullcpp/QhullVertex.h>
//-#include <libqhullcpp/Qhull.h>
#include "interfaces.hh"

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
        } else {
          if (XT::Common::FloatCmp::lt(slope_right_abs, slope_centered_abs))
            ret.set_entry(ii, slope_right.get_entry(ii));
          else
            ret.set_entry(ii, centered_slope.get_entry(ii));
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
  typedef Dune::QuadratureRule<RangeFieldType, dimDomain - 1> QuadratureType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> ResultType;


  template <class... Args>
  explicit LocalCouplingFvOperator(const QuadratureType& quadrature, Args&&... args)
    : quadrature_(quadrature)
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
    for (const auto& quad_point : quadrature_) {
      result += ResultType(numerical_flux_.evaluate(local_functions_tuple_entity,
                                                    local_functions_tuple_neighbor,
                                                    *local_source_entity,
                                                    *local_source_neighbor,
                                                    intersection,
                                                    quad_point.position()))
                * quad_point.weight() * intersection.geometry().integrationElement(quad_point.position());
    }
    local_range_entity.vector().add(result * (1.0 / entity.geometry().volume()));
    local_range_neighbor.vector().add(result * (-1.0 / neighbor.geometry().volume()));
  }

private:
  const QuadratureType& quadrature_;
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
  typedef Dune::QuadratureRule<RangeFieldType, dimDomain - 1> QuadratureType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> ResultType;

  template <class... Args>
  explicit LocalBoundaryFvOperator(const QuadratureType& quadrature, Args&&... args)
    : quadrature_(quadrature)
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
    for (const auto& quad_point : quadrature_) {
      result += ResultType(numerical_flux_.evaluate(
                    local_functions_tuple, *local_source_entity, intersection, quad_point.position()))
                * quad_point.weight() * intersection.geometry().integrationElement(quad_point.position());
    }
    result /= entity.geometry().volume();
    local_range_entity.vector().add(result);
  }

private:
  const QuadratureType& quadrature_;
  const NumericalFluxType numerical_flux_;
}; // class LocalBoundaryFvOperator<...>

template <class SourceType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
class LocalRealizabilityLimiter : public XT::Grid::Functor::Codim0<typename SourceType::SpaceType::GridLayerType>
{
  typedef typename SourceType::SpaceType::GridLayerType GridLayerType;
  typedef typename SourceType::EntityType EntityType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate DomainType;
  typedef typename SourceType::RangeType RangeType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename XT::LA::EigenDenseVector<RangeFieldType> EigenVectorType;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;
  typedef typename std::vector<std::pair<RangeType, RangeFieldType>> PlaneCoefficientsType;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit LocalRealizabilityLimiter(
      const SourceType& source,
      std::vector<std::map<DomainType, RangeType, internal::FieldVectorLess>>& reconstructed_values,
      const BasisFunctionType& basis_functions,
      const QuadratureType& quadrature,
      RangeFieldType epsilon = 1e-14)
    : source_(source)
    , index_set_(source_.space().grid_layer().indexSet())
    , reconstructed_values_(reconstructed_values)
    , basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , epsilon_(epsilon)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    is_instantiated_ = true;
    if (!plane_coefficients_)
      plane_coefficients_ = calculate_plane_coefficients();
  }

  ~LocalRealizabilityLimiter()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_values_[index_set_.index(entity)];

    // get cell average
    const RangeType& u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      // rescale u_l, u_bar
      auto u_l = pair.second;
      auto u_minus_u_bar_l = u_bar - u_l;
      const auto factor = basis_functions_.realizability_limiter_max(u_l, u_bar);
      u_l /= factor;
      u_minus_u_bar_l /= factor;

      for (const auto& coeffs : *plane_coefficients_) {
        const RangeType& a = coeffs.first;
        const RangeFieldType& b = coeffs.second;
        RangeFieldType theta_li = (b - a * u_l) / (a * u_minus_u_bar_l);
        if (XT::Common::FloatCmp::ge(theta_li, -epsilon_) && XT::Common::FloatCmp::le(theta_li, 1.))
          thetas[ll] = std::max(thetas[ll], theta_li);
      } // coeffs
    } // ll
    for (auto& theta : thetas)
      theta = std::min(epsilon_ + theta, 1.);

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    if (XT::Common::FloatCmp::ne(theta_entity, 0.)) {
      for (auto& pair : local_reconstructed_values) {
        auto& u = pair.second;
        auto u_scaled = u;
        u_scaled *= (1 - theta_entity);
        auto u_bar_scaled = u_bar;
        u_bar_scaled *= theta_entity;
        u = u_scaled + u_bar_scaled;
      }
    }
  } // void apply_local(...)

private:
  // calculate half space representation of realizable set
  std::shared_ptr<const PlaneCoefficientsType> calculate_plane_coefficients()
  {
    using orgQhull::Qhull;
    Qhull qhull;
    std::vector<FieldVector<RangeFieldType, dimRange>> points(quadrature_.size() + 1);
    points[0] = FieldVector<RangeFieldType, dimRange>(0);
    size_t ii = 1;
    for (const auto& quad_point : quadrature_)
      points[ii++] = basis_functions_.evaluate(quad_point.position());

    qhull.runQhull("Realizable set", int(dimRange), int(points.size()), &(points[0][0]), "Qt");
    //    qhull.outputQhull("n");
    const auto facet_end = qhull.endFacet();
    std::vector<std::pair<RangeType, RangeFieldType>> plane_coefficients(qhull.facetList().count());
    ii = 0;
    for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ii) {
      for (size_t jj = 0; jj < dimRange; ++jj)
        plane_coefficients[ii].first[jj] = *(facet.hyperplane().coordinates() + jj);
      plane_coefficients[ii].second = -facet.hyperplane().offset();
    }
    return std::make_shared<const PlaneCoefficientsType>(plane_coefficients);
  }

  const SourceType& source_;
  const IndexSetType& index_set_;
  std::vector<std::map<DomainType, RangeType, internal::FieldVectorLess>>& reconstructed_values_;
  const BasisFunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const RangeFieldType epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<const PlaneCoefficientsType> plane_coefficients_;
}; // class LocalRealizabilityLimiter

template <class SourceType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
bool LocalRealizabilityLimiter<SourceType, BasisFunctionType, dimDomain, dimRange>::is_instantiated_ = false;

template <class SourceType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
std::shared_ptr<const typename LocalRealizabilityLimiter<SourceType, BasisFunctionType, dimDomain, dimRange>::
                    PlaneCoefficientsType>
    LocalRealizabilityLimiter<SourceType, BasisFunctionType, dimDomain, dimRange>::plane_coefficients_;

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
class LocalReconstructionFvOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t stencil_size = 2 * polOrder + 1;
  static constexpr std::array<int, 3> stencil = {
      2 * polOrder + 1, dimDomain > 1 ? 2 * polOrder + 1 : 1, dimDomain > 2 ? 2 * polOrder + 1 : 1};
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef typename Dune::XT::Common::FieldVector<RangeFieldType, dimRange> XTFieldVectorType;
  typedef typename XT::LA::EigenDenseVector<RangeFieldType> EigenVectorType;
  typedef typename XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> QuadratureType;
  typedef FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> IntersectionVectorType;
  typedef FieldVector<FieldVector<FieldVector<EigenVectorType, stencil[2]>, stencil[1]>, stencil[0]> ValuesType;
  typedef typename GridLayerType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;

public:
  explicit LocalReconstructionFvOperator(
      const std::vector<EigenVectorType> source_values,
      const AnalyticalFluxType& analytical_flux,
      const BoundaryValueType& boundary_values,
      const GridLayerType& grid_layer,
      const XT::Common::Parameter& param,
      const QuadratureType& quadrature,
      std::vector<std::map<DomainType, RangeType, internal::FieldVectorLess>>& reconstructed_values)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(grid_layer)
    , param_(param)
    , quadrature_(quadrature)
    , reconstructed_values_(reconstructed_values)
  {
    param_.set("boundary", {0.});
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
    const auto& u_entity = XT::LA::internal::FieldVectorToLaVector<EigenVectorType, dimRange>::convert_back(
        values[stencil[0] / 2][stencil[1] / 2][stencil[2] / 2]);
    const FieldVector<FieldMatrix<DomainFieldType, dimRange, dimRange>, dimDomain> jacobians =
        analytical_flux_.local_function(entity)->jacobian_wrt_u(
            entity.geometry().local(entity.geometry().center()), u_entity, param_);

    const auto& entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_values_[entity_index];

    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // get eigenvectors of jacobian
      const auto jacobian =
          XT::LA::internal::FieldMatrixToLaDenseMatrix<EigenMatrixType, dimRange, dimRange>::convert(jacobians[dd]);
      ::Eigen::SelfAdjointEigenSolver<typename EigenMatrixType::BackendType> eigen_solver(jacobian.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvectors_backend = eigen_solver.eigenvectors();
      const auto eigenvectors = EigenMatrixType(eigenvectors_backend.real());
      const auto eigenvectors_inverse = EigenMatrixType(eigenvectors_backend.inverse().real());

      helper<>::reconstruct(
          dd, values, eigenvectors, eigenvectors_inverse, quadrature_, reconstructed_values_map, intersections);
    } // dd

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
                            const EigenMatrixType& eigenvectors,
                            const EigenMatrixType& eigenvectors_inverse,
                            const QuadratureType& /*quadrature*/,
                            std::map<DomainType, RangeType, internal::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<EigenVectorType, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        char_values[ii] = eigenvectors_inverse * values[ii][0][0];

      // reconstruction in x direction
      FieldVector<EigenVectorType, 2> reconstructed_values;
      // quadrature rule containing left and right interface points
      const auto left_and_right_boundary_point = left_right_quadrature();
      slope_reconstruction(char_values, reconstructed_values, left_and_right_boundary_point);

      // convert coordinates on face to local entity coordinates and store
      for (size_t ii = 0; ii < 2; ++ii) {
        // convert back to non-characteristic variables and to FieldVector instead of EigenVector
        const auto value = XT::LA::internal::FieldVectorToLaVector<EigenVectorType, dimRange>::convert_back(
            eigenvectors * reconstructed_values[ii]);
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
                            const EigenMatrixType& eigenvectors,
                            const EigenMatrixType& eigenvectors_inverse,
                            const QuadratureType& quadrature,
                            std::map<DomainType, RangeType, internal::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y-direction accordingly. For that purpose, define new
      // coordinates x', y'. First convert to characteristic variables and reorder x, y to x', y'.
      // Reordering is done such that the indices in char_values are in the order y', x'.
      FieldVector<FieldVector<EigenVectorType, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          if (dd == 0)
            char_values[jj][ii] = eigenvectors_inverse * values[ii][jj][0];
          else if (dd == 1)
            char_values[ii][jj] = eigenvectors_inverse * values[ii][jj][0];
        }
      }

      // reconstruction in x' direction
      // first index: dimension 2 for left/right interface
      // second index: y' direction
      FieldVector<FieldVector<EigenVectorType, stencil_size>, 2> x_reconstructed_values;
      std::vector<EigenVectorType> result(2);
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
      FieldVector<std::vector<EigenVectorType>, 2> reconstructed_values(
          (std::vector<EigenVectorType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        slope_reconstruction(x_reconstructed_values[ii], reconstructed_values[ii], quadrature);

      // convert coordinates on face to local entity coordinates and store
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          // convert back to non-characteristic variables and to FieldVector instead of EigenVector
          const auto value = XT::LA::internal::FieldVectorToLaVector<EigenVectorType, dimRange>::convert_back(
              eigenvectors * reconstructed_values[ii][jj]);
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
                            const EigenMatrixType& eigenvectors,
                            const EigenMatrixType& eigenvectors_inverse,
                            const QuadratureType& quadrature,
                            std::map<DomainType, RangeType, internal::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y- and z-direction accordingly. For that purpose, define new
      // coordinates x', y', z'. First convert to characteristic variables and reorder x, y, z to x', y', z'.
      // Reordering is done such that the indices in char_values are in the order z', y', x'
      FieldVector<FieldVector<FieldVector<EigenVectorType, stencil_size>, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            if (dd == 0)
              char_values[kk][jj][ii] = eigenvectors_inverse * values[ii][jj][kk];
            else if (dd == 1)
              char_values[ii][kk][jj] = eigenvectors_inverse * values[ii][jj][kk];
            else if (dd == 2)
              char_values[jj][ii][kk] = eigenvectors_inverse * values[ii][jj][kk];
          }
        }
      }

      // reconstruction in x' direction
      // first index: z' direction
      // second index: dimension 2 for left/right interface
      // third index: y' direction
      FieldVector<FieldVector<FieldVector<EigenVectorType, stencil_size>, 2>, stencil_size> x_reconstructed_values;
      std::vector<EigenVectorType> result(2);
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
      FieldVector<std::vector<FieldVector<EigenVectorType, stencil_size>>, 2> y_reconstructed_values(
          (std::vector<FieldVector<EigenVectorType, stencil_size>>(num_quad_points)));
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
      FieldVector<std::vector<std::vector<EigenVectorType>>, 2> reconstructed_values(
          std::vector<std::vector<EigenVectorType>>(num_quad_points, std::vector<EigenVectorType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < num_quad_points; ++jj)
          slope_reconstruction(y_reconstructed_values[ii][jj], reconstructed_values[ii][jj], quadrature);

      // convert coordinates on face to local entity coordinates and store
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < num_quad_points; ++kk) {
            // convert back to non-characteristic variables and to FieldVector instead of EigenVector
            const auto value = XT::LA::internal::FieldVectorToLaVector<EigenVectorType, dimRange>::convert_back(
                eigenvectors * reconstructed_values[ii][jj][kk]);
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

    static void apply(const std::vector<EigenVectorType>& source_values,
                      const BoundaryValueType& boundary_values,
                      FieldVector<FieldVector<FieldVector<EigenVectorType, stencil_z>, stencil_y>, stencil_x>& values,
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
                  XT::LA::internal::FieldVectorToLaVector<EigenVectorType, dimRange>::convert(boundary_value);
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
      return !(std::abs(offsets[dir / 2]) < stencil[dir / 2] / 2);
    }
  }; // class StencilIterator<...>

  static void slope_reconstruction(const FieldVector<EigenVectorType, stencil_size>& cell_values,
                                   std::vector<EigenVectorType>& result,
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
    const auto slope_centered = (u_right - u_left) * 0.5;
    const auto slope =
        internal::ChooseLimiter<slope_limiter, EigenVectorType>::limit(slope_left, slope_right, slope_centered);
    for (size_t ii = 0; ii < points.size(); ++ii)
      result[ii] = u_entity + slope * (points[ii] - 0.5);
  }

  static void slope_reconstruction(const FieldVector<EigenVectorType, stencil_size>& cell_values,
                                   FieldVector<EigenVectorType, 2>& result,
                                   const QuadratureType& quadrature)
  {
    std::vector<EigenVectorType> result_vec(2);
    slope_reconstruction(cell_values, result_vec, quadrature);
    std::copy(result_vec.begin(), result_vec.end(), result.begin());
  }


  const std::vector<EigenVectorType> source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const QuadratureType quadrature_;
  std::vector<std::map<DomainType, RangeType, internal::FieldVectorLess>>& reconstructed_values_;
}; // class LocalReconstructionFvOperator

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
constexpr std::array<int, 3>
    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
        stencil;

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_FV_HH
