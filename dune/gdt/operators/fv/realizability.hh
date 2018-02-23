// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
#define DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker/functors.hh>

#if HAVE_QHULL
#include <dune/xt/common/disable_warnings.hh>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <dune/xt/common/reenable_warnings.hh>
#endif // HAVE_QHULL

#if HAVE_LPSOLVE
#include <dune/xt/common/disable_warnings.hh>
namespace lpsolve {
#include <lpsolve/lp_lib.h>
}
#include <dune/xt/common/reenable_warnings.hh>
#endif // HAVE_LPSOLVE

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType, class BasisfunctionType>
class RealizabilityLimiterBase
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
public:
  typedef typename DiscreteFunctionType::SpaceType::GridLayerType GridLayerType;
  typedef typename DiscreteFunctionType::EntityType EntityType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate DomainType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  static const size_t dimDomain = DiscreteFunctionType::dimDomain;
  static const size_t dimRange = DiscreteFunctionType::dimRange;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;

  explicit RealizabilityLimiterBase(const BasisfunctionType& basis_functions,
                                    const QuadratureType& quadrature,
                                    RangeFieldType epsilon)
    : basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , num_quad_points_(quadrature_.size())
    , epsilon_(epsilon)
    , M_(quadrature_.size())
  {
    for (size_t ii = 0; ii < quadrature_.size(); ++ii)
      M_[ii] = basis_functions_.evaluate(quadrature_[ii].position());
  }

  void set_index_set(const IndexSetType* index_set)
  {
    index_set_ = index_set;
  }

  void set_source_values(const std::vector<RangeType>* source_values)
  {
    source_values_ = source_values;
  }

  void set_reconstructed_values(
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>* reconstructed_values)
  {
    reconstructed_values_ = reconstructed_values;
  }

protected:
  const BasisfunctionType& basis_functions_;
  const std::vector<RangeType>* source_values_;
  const IndexSetType* index_set_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>* reconstructed_values_;
  const QuadratureType& quadrature_;
  const size_t num_quad_points_;
  const RangeFieldType epsilon_;
  std::vector<RangeType> M_;
};

template <class DiscreteFunctionType, class BasisfunctionType = int>
class NonLimitingRealizabilityLimiter : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::QuadratureType;
  using typename BaseType::RangeFieldType;

  explicit NonLimitingRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                           const QuadratureType& quadrature,
                                           RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& /*entity*/)
  {
  }
}; // class NonLimitingRealizabilityLimiter

template <class DiscreteFunctionType, class BasisfunctionType>
class RelativePositivityLocalRealizabilityLimiter
    : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit RelativePositivityLocalRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                                       const QuadratureType& quadrature,
                                                       RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& entity)
  {
    const auto entity_index = index_set_->index(entity);
    auto& local_reconstructed_values = (*reconstructed_values_)[entity_index];

    // get cell average
    const RangeType& u_bar = (*source_values_)[entity_index];

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), 0.);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      for (size_t ii = 0; ii < u_l.size(); ++ii)
        if (XT::Common::FloatCmp::ne(u_l[ii], u_bar[ii])) {
          auto theta_ii = u_l[ii] / (u_l[ii] - u_bar[ii]) + epsilon_;
          if (XT::Common::FloatCmp::le(theta_ii, 1.))
            thetas[ll] = std::max(thetas[ll], theta_ii);
        }
    } // ll
    for (auto& theta : thetas)
      theta = std::min(theta, 1.);

    auto theta_entity_it = std::max_element(thetas.begin(), thetas.end());
    auto theta_entity = theta_entity_it == thetas.end() ? 0. : *theta_entity_it;
    if (theta_entity > 0.) {
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
  using BaseType::source_values_;
  using BaseType::index_set_;
  using BaseType::reconstructed_values_;
  using BaseType::epsilon_;
}; // class RelativePositivityLocalRealizabilityLimiter

template <class DiscreteFunctionType, class BasisfunctionType>
class PositivityLocalRealizabilityLimiter : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit PositivityLocalRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                               const QuadratureType& quadrature,
                                               RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& entity)
  {
    const auto entity_index = index_set_->index(entity);
    auto& local_reconstructed_values = (*reconstructed_values_)[entity_index];

    // get cell average
    const RangeType& u_bar = (*source_values_)[entity_index];

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), 0.);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      for (size_t ii = 0; ii < u_l.size(); ++ii)
        if (u_l[ii] < epsilon_ && u_bar[ii] > epsilon_)
          thetas[ll] = std::max(thetas[ll], (epsilon_ - u_l[ii]) / (u_bar[ii] - u_l[ii]));
    } // ll
    for (auto& theta : thetas)
      theta = std::min(theta, 1.);

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    if (theta_entity > 0.) {
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
  using BaseType::source_values_;
  using BaseType::index_set_;
  using BaseType::reconstructed_values_;
  using BaseType::epsilon_;
}; // class PositivityLocalRealizabilityLimiter

template <class DiscreteFunctionType, class BasisfunctionType>
class DgLocalRealizabilityLimiter : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit DgLocalRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                       const QuadratureType& quadrature,
                                       RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
    , triangulation_(basis_functions.triangulation())
  {
  }

  void apply_local(const EntityType& entity)
  {
    const auto entity_index = index_set_->index(entity);
    auto& local_reconstructed_values = (*reconstructed_values_)[entity_index];

    // get cell average
    const RangeType u_bar = (*source_values_)[entity_index];

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      FieldVector<RangeFieldType, dimRange / 2> local_thetas(-epsilon_);
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        const auto& u0 = u_l[2 * ii];
        const auto& u1 = u_l[2 * ii + 1];
        const auto& ubar0 = u_bar[2 * ii];
        const auto& ubar1 = u_bar[2 * ii + 1];
        const auto& vj = triangulation_[ii];
        const auto& vjplus1 = triangulation_[ii];
        FieldVector<RangeFieldType, 3> thetas_ii;
        thetas_ii[0] = u0 / (ubar0 - u0);
        thetas_ii[1] = (u0 * vj - u1) / ((ubar1 - u1) - (ubar0 - u0) * vj);
        thetas_ii[2] = (u0 * vjplus1 - u1) / ((ubar1 - u1) - (ubar0 - u0) * vjplus1);
        for (size_t kk = 0; kk < 3; ++kk) {
          if (thetas_ii[kk] < -epsilon_ || thetas_ii[kk] > 1.)
            thetas_ii[kk] = -epsilon_;
          local_thetas[ii] = std::max(local_thetas[ii], thetas_ii[kk]);
        } // kk
      } // ii
      thetas[ll] = *std::max_element(local_thetas.begin(), local_thetas.end());
    } // ll
    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    theta_entity = std::min(theta_entity + epsilon_, 1.);
    if (theta_entity > 0.) {
      std::cout << "Realizability limiter applied, theta = " << theta_entity << std::endl;
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
  using BaseType::source_values_;
  using BaseType::index_set_;
  using BaseType::reconstructed_values_;
  using BaseType::epsilon_;
  typename BasisfunctionType::TriangulationType triangulation_;
}; // class DgLocalRealizabilityLimiter

#if HAVE_QHULL

template <class DiscreteFunctionType, class BasisfunctionType>
class ConvexHullLocalRealizabilityLimiter : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  typedef typename std::vector<std::pair<RangeType, RangeFieldType>> PlaneCoefficientsType;

  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit ConvexHullLocalRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                               const QuadratureType& quadrature,
                                               RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    is_instantiated_ = true;
    if (!plane_coefficients_)
      calculate_plane_coefficients();
  }

  ~ConvexHullLocalRealizabilityLimiter()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    auto entity_index = index_set_->index(entity);
    auto& local_reconstructed_values = (*reconstructed_values_)[entity_index];

    // get cell average
    const RangeType& u_bar = (*source_values_)[entity_index];

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      // rescale u_l, u_bar
      auto u_l = pair.second;
      auto u_bar_minus_u_l = u_bar - u_l;
      const auto factor = basis_functions_.realizability_limiter_max(u_l, u_bar);
      u_l /= factor;
      u_bar_minus_u_l /= factor;

      for (const auto& coeffs : *plane_coefficients_) {
        const RangeType& a = coeffs.first;
        const RangeFieldType& b = coeffs.second;
        RangeFieldType theta_li = (b - a * u_l) / (a * u_bar_minus_u_l);
        if (XT::Common::FloatCmp::le(theta_li, 1.))
          thetas[ll] = std::max(thetas[ll], theta_li);
      } // coeffs
    } // ll
    for (auto& theta : thetas)
      theta = std::min(epsilon_ + theta, 1.);

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    if (theta_entity > 0.) {
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
  void calculate_plane_coefficients()
  {
    using orgQhull::Qhull;
    Qhull qhull;
    std::vector<FieldVector<RangeFieldType, dimRange>> points(quadrature_.size() + 1);
    points[0] = FieldVector<RangeFieldType, dimRange>(0);
    size_t ii = 1;
    for (const auto& quad_point : quadrature_)
      points[ii++] = basis_functions_.evaluate(quad_point.position());

    std::cout << "Starting qhull..." << std::endl;
    qhull.runQhull("Realizable set", int(dimRange), int(points.size()), &(points[0][0]), "Qt T1");
    std::cout << "qhull done" << std::endl;
    //    qhull.outputQhull("n");
    const auto facet_end = qhull.endFacet();
    plane_coefficients_ = std::make_shared<PlaneCoefficientsType>(qhull.facetList().count());
    ii = 0;
    for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ii) {
      for (size_t jj = 0; jj < dimRange; ++jj)
        (*plane_coefficients_)[ii].first[jj] = *(facet.hyperplane().coordinates() + jj);
      (*plane_coefficients_)[ii].second = -facet.hyperplane().offset();
    }
  }

  using BaseType::basis_functions_;
  using BaseType::source_values_;
  using BaseType::index_set_;
  using BaseType::reconstructed_values_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<PlaneCoefficientsType> plane_coefficients_;
}; // class ConvexHullLocalRealizabilityLimiter

template <class DiscreteFunctionType, class BasisfunctionType>
bool ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType>::is_instantiated_ = false;

template <class DiscreteFunctionType, class BasisfunctionType>
std::shared_ptr<
    typename ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType>::PlaneCoefficientsType>
    ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType>::plane_coefficients_;

template <class DiscreteFunctionType, class BasisfunctionType>
class DgConvexHullLocalRealizabilityLimiter : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t block_size = (dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  typedef FieldVector<RangeFieldType, block_size> BlockRangeType;
  typedef typename std::vector<std::pair<BlockRangeType, RangeFieldType>> BlockPlaneCoefficientsType;
  typedef FieldVector<BlockPlaneCoefficientsType, num_blocks> PlaneCoefficientsType;

  explicit DgConvexHullLocalRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                                 const QuadratureType& quadrature,
                                                 RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    is_instantiated_ = true;
    if (!plane_coefficients_)
      calculate_plane_coefficients();
  }

public:
  ~DgConvexHullLocalRealizabilityLimiter()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    const auto entity_index = index_set_->index(entity);
    auto& local_reconstructed_values = (*reconstructed_values_)[entity_index];
    const RangeType& u_bar = (*source_values_)[entity_index];

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      // rescale u_l, u_bar
      auto u_l = pair.second;
      if (XT::Common::FloatCmp::eq(u_l, u_bar))
        continue;
      auto u_bar_minus_u_l = u_bar - u_l;
      const auto factor = basis_functions_.realizability_limiter_max(u_l, u_bar);
      u_l /= factor;
      u_bar_minus_u_l /= factor;

      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t offset = jj * block_size;
        for (const auto& coeffs : (*plane_coefficients_)[jj]) {
          const BlockRangeType& a = coeffs.first;
          const RangeFieldType& b = coeffs.second;
          RangeFieldType a_dot_u_l_block(0);
          RangeFieldType a_dot_u_bar_minus_u_l_block(0);
          for (size_t mm = 0; mm < block_size; ++mm) {
            const size_t index = offset + mm;
            a_dot_u_l_block += u_l[index] * a[mm];
            a_dot_u_bar_minus_u_l_block += u_bar_minus_u_l[index] * a[mm];
          }
          RangeFieldType theta_li = (b - a_dot_u_l_block) / a_dot_u_bar_minus_u_l_block;
          if (XT::Common::FloatCmp::le(theta_li, 1.))
            thetas[ll] = std::max(thetas[ll], theta_li);
        } // coeffs
      } // jj
    } // ll
    for (auto& theta : thetas)
      theta = std::min(epsilon_ + theta, 1.);

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    if (theta_entity > 0.) {
      std::cout << "realizability limited with theta = " << theta_entity << std::endl;
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
  void calculate_plane_coefficients()
  {
    plane_coefficients_ = std::make_shared<PlaneCoefficientsType>();
    using orgQhull::Qhull;
    FieldVector<std::vector<FieldVector<RangeFieldType, block_size>>, num_blocks> points;
    FieldVector<QuadratureType, num_blocks> blocked_quadrature;
    for (const auto& quad_point : quadrature_) {
      const auto face_indices = basis_functions_.get_face_indices(quad_point.position());
      const size_t num_adjacent_faces = face_indices.size();
      for (const auto& kk : face_indices)
        blocked_quadrature[kk].emplace_back(quad_point.position(), quad_point.weight() / num_adjacent_faces);
    } // ii
    size_t num_faces;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      points[jj].resize(blocked_quadrature[jj].size() + 1);
      for (size_t ii = 0; ii < blocked_quadrature[jj].size(); ++ii) {
        const auto val = basis_functions_.evaluate(blocked_quadrature[jj][ii].position(), false, num_faces);
        for (size_t ll = 0; ll < block_size; ++ll)
          points[jj][ii][ll] = val[block_size * jj + ll];
      } // ii
      points[jj][blocked_quadrature[jj].size()] = FieldVector<RangeFieldType, block_size>(0.);
      Qhull qhull;
      std::cout << "Starting qhull..." << std::endl;
      qhull.runQhull("Realizable set", int(block_size), int(points[jj].size()), &(points[jj][0][0]), "Qt T1");
      std::cout << "qhull done" << std::endl;
      const auto facet_end = qhull.endFacet();
      BlockPlaneCoefficientsType block_plane_coefficients(qhull.facetList().count());
      size_t ii = 0;
      for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ii) {
        for (size_t ll = 0; ll < block_size; ++ll)
          block_plane_coefficients[ii].first[ll] = *(facet.hyperplane().coordinates() + ll);
        block_plane_coefficients[ii].second = -facet.hyperplane().offset();
      }
      (*plane_coefficients_)[jj] = block_plane_coefficients;
    } // jj
  }

  using BaseType::basis_functions_;
  using BaseType::source_values_;
  using BaseType::index_set_;
  using BaseType::reconstructed_values_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<PlaneCoefficientsType> plane_coefficients_;
}; // class ConvexHullLocalRealizabilityLimiter

template <class DiscreteFunctionType, class BasisfunctionType>
bool DgConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType>::is_instantiated_ = false;

template <class DiscreteFunctionType, class BasisfunctionType>
std::shared_ptr<
    typename DgConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType>::PlaneCoefficientsType>
    DgConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType>::plane_coefficients_;

#else // HAVE_QHULL

template <class DiscreteFunctionType, class BasisfunctionType>
class ConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionType>::value, "You are missing Qhull!");
};

template <class DiscreteFunctionType, class BasisfunctionType>
class DgConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionType>::value, "You are missing Qhull!");
};

#endif // HAVE_QHULL

#if HAVE_LPSOLVE

template <class DiscreteFunctionType, class BasisfunctionType>
class LPLocalRealizabilityLimiter : public RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType>
{
  typedef RealizabilityLimiterBase<DiscreteFunctionType, BasisfunctionType> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  explicit LPLocalRealizabilityLimiter(const BasisfunctionType& basis_functions,
                                       const QuadratureType& quadrature,
                                       RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& entity)
  {
    const auto entity_index = index_set_->index(entity);
    auto& local_reconstructed_values = (*reconstructed_values_)[entity_index];

    // get cell average
    const RangeType& u_bar = (*source_values_)[entity_index];

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      // if (XT::Common::FloatCmp::eq(u_bar, u_l) || basis_functions_.obviously_realizable(u_l)) {
      if (XT::Common::FloatCmp::eq(u_bar, u_l)) {
        thetas[ll] = 0.;
        continue;
      }

      // solve LP:
      // min \theta s.t.
      // (\sum x_i v_i) + \theta (u - \bar{u}) = u
      // x_i, \theta >= 0
      auto theta = solve_linear_program(u_bar, u_l);
      theta = theta + epsilon_;
      if (theta > 1.)
        theta = 1.;
      thetas[ll] = theta;
    } // ll

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    if (theta_entity > 0.)
      limit_realizability(local_reconstructed_values, u_bar, theta_entity);
  } // void apply_local(...)

private:
  void limit_realizability(std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& values,
                           const RangeType& u_bar,
                           const RangeFieldType theta)
  {
    for (auto& pair : values) {
      auto& u = pair.second;
      u = convex_combination(u, u_bar, theta);
    }
  }

  RangeType convex_combination(const RangeType& u, const RangeType& u_bar, const RangeFieldType& theta)
  {
    if (XT::Common::FloatCmp::eq(theta, 0.))
      return u;
    else {
      RangeType u_scaled = u;
      u_scaled *= 1 - theta;
      RangeType u_bar_scaled = u_bar;
      u_bar_scaled *= theta;
      return u_scaled + u_bar_scaled;
    }
  }

  //  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l, const size_t index)
  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l)
  {
    //        auto& regularization_parameter = regularization_parameters_[index];
    //        auto& regularization_time = regularization_times_[index];
    //        if (XT::Common::FloatCmp::ne(regularization_time, t)) {
    //            regularization_time_ = t;
    //            regularization_parameter_ = 0.;
    //        }
    RangeFieldType theta;

    //    const auto r_max = r_sequence_.back();
    //    for (const auto& r : r_sequence_) {
    //        regularization_parameter = std::max(regularization_parameter, r);

    typename lpsolve::lprec* lp;

    const auto u_l_minus_u_bar = u_l - u_bar;

    // We start with creating a model with dimRange+1 rows and num_quad_points+1 columns */
    constexpr int num_rows = int(dimRange);
    int num_cols = int(num_quad_points_ + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
    lp = lpsolve::make_lp(num_rows, num_cols);
    if (!lp)
      DUNE_THROW(Dune::MathError, "Couldn't construct linear program");

    /* let us name our variables. Not required, but can be useful for debugging */
    for (int ii = 1; ii <= num_cols - 1; ++ii) {
      std::string ii_string = XT::Common::to_string(ii);
      std::string name_string = "x";
      for (size_t jj = 5 - ii_string.size(); jj > 0; --jj)
        name_string += "0";
      name_string += ii_string;
      lpsolve::set_col_name(lp, ii, &name_string[0]);
    }

    std::string name_string = "theta ";
    lpsolve::set_col_name(lp, num_cols, &name_string[0]);

    // In the call to set_column, the first entry (row 0) is the value of the objective function
    // (c_i in the objective function c^T x), the other entries are the entries of the i-th column
    // in the constraints matrix. The 0-th column is the rhs vector. The entry (0, 0) corresponds
    // to the initial value of the objective function.
    std::array<REAL, num_rows + 1> column;

    // set rhs (column 0)
    column[0] = 0.;
    std::copy(u_l.begin(), u_l.end(), column.begin() + 1);
    lpsolve::set_rh_vec(lp, column.data());
    lpsolve::set_rh(lp, 0, column[0]);

    // set columns for quadrature points
    column[0] = 0.;
    for (int ii = 0; ii < int(M_.size()); ++ii) {
      const auto& v_i = M_[ii];
      std::copy(v_i.begin(), v_i.end(), column.begin() + 1);
      lpsolve::set_column(lp, ii + 1, column.data());
    }

    // set theta column
    column[0] = 1.;
    std::copy(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end(), column.begin() + 1);
    lpsolve::set_column(lp, num_cols, column.data());

    // set all contraints to equality constraints
    for (int ii = 1; ii <= num_rows; ++ii)
      lpsolve::set_constr_type(lp, ii, EQ);

    // Set bounds for variables. 0 <= x <= inf is the default for all variables.
    // The bounds for theta should be [0,1], but we allow allow theta to be -0.1
    // so we can check if u_l is on the boundary and if so, move it a little away
    // from the boundary
    lpsolve::set_bounds(lp, num_cols, -0.1, 1.);

    default_basis(lp);

    /* I only want to see important messages on screen while solving */
    lpsolve::set_verbose(lp, IMPORTANT);

    /* Now let lpsolve calculate a solution */
    const auto solve_status = lpsolve::solve(lp);
    theta = lpsolve::get_objective(lp);
    if (solve_status != OPTIMAL) {
      lpsolve::write_LP(lp, stdout);
      std::cout << solve_status << std::endl;
      //            if (r == r_max)
      DUNE_THROW(Dune::MathError, "An unexpected error occured while solving the linear program");
      //            break;
    }

    /* free allocated memory */
    if (lp)
      lpsolve::delete_lp(lp);
    //        }

    return theta;
  }

  using BaseType::source_values_;
  using BaseType::index_set_;
  using BaseType::reconstructed_values_;
  using BaseType::num_quad_points_;
  using BaseType::epsilon_;
  using BaseType::M_;
}; // class LPLocalRealizabilityLimiter

#else // HAVE_LPSOLVE

template <class DiscreteFunctionType, class BasisfunctionType>
class LPLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionType>::value, "You are missing LPSolve!");
};

#endif // HAVE_LPSOLVE


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
