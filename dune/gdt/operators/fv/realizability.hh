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

#include <dune/xt/grid/walker.hh>

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

#include <dune/gdt/operators/interfaces.hh>

#include "reconstructed_function.hh"

namespace Dune {
namespace GDT {


template <class DiscreteFunctionImp, class BasisfunctionImp>
class LocalRealizabilityLimiterBase
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionImp::SpaceType::GridLayerType>
{
public:
  using DiscreteFunctionType = DiscreteFunctionImp;
  using BasisfunctionType = BasisfunctionImp;
  typedef typename DiscreteFunctionType::SpaceType::GridLayerType GridLayerType;
  typedef typename DiscreteFunctionType::EntityType EntityType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate DomainType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = DiscreteFunctionType::dimDomain;
  static const size_t dimRange = DiscreteFunctionType::dimRange;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;
  using ReconstructedFunctionType =
      ReconstructedLocalizableFunction<GridLayerType, DomainFieldType, dimDomain, RangeFieldType, dimRange>;

  explicit LocalRealizabilityLimiterBase(const DiscreteFunctionType& source,
                                         ReconstructedFunctionType& reconstructed_function,
                                         const BasisfunctionType& basis_functions,
                                         const QuadratureType& quadrature,
                                         const RangeFieldType epsilon,
                                         const std::vector<RangeType>& basis_values)
    : source_(source)
    , reconstructed_function_(reconstructed_function)
    , basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , epsilon_(epsilon)
    , basis_values_(basis_values)
  {
  }

protected:
  const DiscreteFunctionType& source_;
  ReconstructedFunctionType& reconstructed_function_;
  const BasisfunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const RangeFieldType epsilon_;
  const std::vector<RangeType>& basis_values_;
};

template <class DiscreteFunctionImp, class BasisfunctionImp = int>
class NonLimitingLocalRealizabilityLimiter : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::QuadratureType;
  using typename BaseType::RangeFieldType;

  explicit NonLimitingLocalRealizabilityLimiter(const BasisfunctionImp& basis_functions,
                                                const QuadratureType& quadrature,
                                                RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& /*entity*/)
  {
  }
}; // class NonLimitingLocalRealizabilityLimiter

template <class DiscreteFunctionImp, class BasisfunctionImp>
class RelativePositivityLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

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
  explicit RelativePositivityLocalRealizabilityLimiter(const BasisfunctionImp& basis_functions,
                                                       const QuadratureType& quadrature,
                                                       RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

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
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::epsilon_;
}; // class RelativePositivityLocalRealizabilityLimiter

template <class DiscreteFunctionImp, class BasisfunctionImp>
class PositivityLocalRealizabilityLimiter : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

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
  explicit PositivityLocalRealizabilityLimiter(const BasisfunctionImp& basis_functions,
                                               const QuadratureType& quadrature,
                                               RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

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
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::epsilon_;
}; // class PositivityLocalRealizabilityLimiter

template <class DiscreteFunctionImp, class BasisfunctionImp>
class DgLocalRealizabilityLimiter : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::EntityType;
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::ReconstructedFunctionType;
  using typename BaseType::BasisfunctionType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  explicit DgLocalRealizabilityLimiter(const DiscreteFunctionType& source,
                                       ReconstructedFunctionType& reconstructed_function,
                                       const BasisfunctionType& basis_functions,
                                       const QuadratureType& quadrature,
                                       const RangeFieldType epsilon,
                                       const std::vector<RangeType>& basis_values)
    : BaseType(source, reconstructed_function, basis_functions, quadrature, epsilon, basis_values)
    , triangulation_(basis_functions.triangulation())
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

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
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::epsilon_;
  typename BasisfunctionImp::TriangulationType triangulation_;
}; // class DgLocalRealizabilityLimiter

#if HAVE_QHULL

template <class DiscreteFunctionImp, class BasisfunctionImp>
class ConvexHullLocalRealizabilityLimiter : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

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
  explicit ConvexHullLocalRealizabilityLimiter(const BasisfunctionImp& basis_functions,
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
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

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
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<PlaneCoefficientsType> plane_coefficients_;
}; // class ConvexHullLocalRealizabilityLimiter

template <class DiscreteFunctionImp, class BasisfunctionImp>
bool ConvexHullLocalRealizabilityLimiter<DiscreteFunctionImp, BasisfunctionImp>::is_instantiated_ = false;

template <class DiscreteFunctionImp, class BasisfunctionImp>
std::shared_ptr<
    typename ConvexHullLocalRealizabilityLimiter<DiscreteFunctionImp, BasisfunctionImp>::PlaneCoefficientsType>
    ConvexHullLocalRealizabilityLimiter<DiscreteFunctionImp, BasisfunctionImp>::plane_coefficients_;

template <class DiscreteFunctionImp, class BasisfunctionImp>
class DgConvexHullLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

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

  explicit DgConvexHullLocalRealizabilityLimiter(const BasisfunctionImp& basis_functions,
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
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

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
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<PlaneCoefficientsType> plane_coefficients_;
}; // class ConvexHullLocalRealizabilityLimiter

template <class DiscreteFunctionImp, class BasisfunctionImp>
bool DgConvexHullLocalRealizabilityLimiter<DiscreteFunctionImp, BasisfunctionImp>::is_instantiated_ = false;

template <class DiscreteFunctionImp, class BasisfunctionImp>
std::shared_ptr<
    typename DgConvexHullLocalRealizabilityLimiter<DiscreteFunctionImp, BasisfunctionImp>::PlaneCoefficientsType>
    DgConvexHullLocalRealizabilityLimiter<DiscreteFunctionImp, BasisfunctionImp>::plane_coefficients_;

#else // HAVE_QHULL

template <class DiscreteFunctionImp, class BasisfunctionImp>
class ConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionImp>::value, "You are missing Qhull!");
};

template <class DiscreteFunctionImp, class BasisfunctionImp>
class DgConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionImp>::value, "You are missing Qhull!");
};

#endif // HAVE_QHULL

#if HAVE_LPSOLVE

template <class DiscreteFunctionImp, class BasisfunctionImp>
class LPLocalRealizabilityLimiter : public LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp>
{
  typedef LocalRealizabilityLimiterBase<DiscreteFunctionImp, BasisfunctionImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  explicit LPLocalRealizabilityLimiter(const BasisfunctionImp& basis_functions,
                                       const QuadratureType& quadrature,
                                       RangeFieldType epsilon)
    : BaseType(basis_functions, quadrature, epsilon)
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

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
    int num_cols = int(quadrature_.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
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
    for (int ii = 0; ii < int(basis_values_.size()); ++ii) {
      const auto& v_i = basis_values_[ii];
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

  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  using BaseType::basis_values_;
}; // class LPLocalRealizabilityLimiter

#else // HAVE_LPSOLVE

template <class DiscreteFunctionImp, class BasisfunctionImp>
class LPLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionImp>::value, "You are missing LPSolve!");
};

#endif // HAVE_LPSOLVE


template <class LocalRealizabilityLimiterImp, class Traits>
class RealizabilityLimiter;


namespace internal {


template <class LocalRealizabilityLimiterImp>
struct RealizabilityLimiterTraits
{
  using LocalRealizabilityLimiterType = LocalRealizabilityLimiterImp;
  using BasisfunctionType = typename LocalRealizabilityLimiterImp::BasisfunctionType;
  using QuadratureType = typename LocalRealizabilityLimiterImp::QuadratureType;
  using RangeFieldType = typename LocalRealizabilityLimiterImp::RangeFieldType;
  using FieldType = RangeFieldType;
  using JacobianType = NoJacobian;
  using RangeType = typename LocalRealizabilityLimiterImp::RangeType;
  using ReconstructedFunctionType = typename LocalRealizabilityLimiterImp::ReconstructedFunctionType;
  using derived_type = RealizabilityLimiter<LocalRealizabilityLimiterType, RealizabilityLimiterTraits>;
};


} // namespace internal


template <class LocalRealizabilityLimiterImp,
          class Traits = internal::RealizabilityLimiterTraits<LocalRealizabilityLimiterImp>>
class RealizabilityLimiter : public OperatorInterface<Traits>
{
public:
  using LocalRealizabilityLimiterType = typename Traits::LocalRealizabilityLimiterType;
  using BasisfunctionType = typename Traits::BasisfunctionType;
  using QuadratureType = typename Traits::QuadratureType;
  using RangeFieldType = typename Traits::RangeFieldType;
  using RangeType = typename Traits::RangeType;
  using ReconstructedFunctionType = typename Traits::ReconstructedFunctionType;

  RealizabilityLimiter(const BasisfunctionType& basis_functions,
                       const QuadratureType& quadrature,
                       const RangeFieldType epsilon)
    : basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , epsilon_(epsilon)
    , basis_values_(quadrature_.size())
  {
    for (size_t ii = 0; ii < quadrature_.size(); ++ii)
      basis_values_[ii] = basis_functions_.evaluate(quadrature_[ii].position());
  }

  template <class SourceType>
  void apply(const SourceType& source, ReconstructedFunctionType& range, const XT::Common::Parameter& /*param*/) const
  {
    static_assert(is_discrete_function<SourceType>::value,
                  "SourceType has to be derived from DiscreteFunction (use the non-reconstructed values!)");
    LocalRealizabilityLimiterType local_realizability_limiter(
        source, range, basis_functions_, quadrature_, epsilon_, basis_values_);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(source.space().grid_layer());
    walker.append(local_realizability_limiter);
    walker.walk(true);
  } // void apply(...)

private:
  const BasisfunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const RangeFieldType epsilon_;
  std::vector<RangeType> basis_values_;
}; // class RealizabilityLimiter<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
