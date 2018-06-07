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

#include "config.h"

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/walker.hh>

#if HAVE_QHULL
#include <dune/xt/common/disable_warnings.hh>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <dune/xt/common/reenable_warnings.hh>
#endif // HAVE_QHULL

#if HAVE_CLP
#include <coin/ClpSimplex.hpp>
#endif // HAVE_CLP

#include <dune/xt/common/lpsolve.hh>

#include <dune/gdt/operators/fv/reconstruction/reconstructed_function.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/local/fluxes/entropybased.hh>

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class LocalRealizabilityLimiterBase
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionImp::SpaceType::GridLayerType>
{
public:
  using AnalyticalFluxType = AnalyticalFluxImp;
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
  using EntropyFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, DiscreteFunctionType>;

  LocalRealizabilityLimiterBase(const AnalyticalFluxType& analytical_flux,
                                const DiscreteFunctionType& source,
                                ReconstructedFunctionType& reconstructed_function,
                                const BasisfunctionType& basis_functions,
                                const QuadratureType& quadrature,
                                const RangeFieldType epsilon,
                                const std::vector<RangeType>& basis_values,
                                const XT::Common::Parameter& param)
    : analytical_flux_(analytical_flux)
    , source_(source)
    , reconstructed_function_(reconstructed_function)
    , basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , epsilon_(epsilon)
    , basis_values_(basis_values)
    , param_(param)
  {
    param_.set("boundary", {0.});
  }

  virtual ~LocalRealizabilityLimiterBase()
  {
  }

  void apply_limiter(const typename AnalyticalFluxType::EntityType& entity,
                     const RangeFieldType theta_entity,
                     std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& local_reconstructed_values,
                     const RangeType& u_bar)
  {
    assert(dynamic_cast<const EntropyFluxType*>(&analytical_flux_) != nullptr
           && "analytical_flux_ has to be derived from EntropyBasedLocalFlux");
    if (theta_entity + epsilon_ > 0.) {
      std::cout << "limited with theta: " << theta_entity << " and epsilon " << epsilon_ << std::endl;
      auto theta = theta_entity + epsilon_;
      if (theta > 1.)
        theta = 1.;
      for (auto& pair : local_reconstructed_values) {
        auto& u = pair.second;
        u = convex_combination(u, u_bar, theta);
      }
    }
    // solve optimization problem
    for (auto& pair : local_reconstructed_values) {
      const auto x_in_inside_coords = entity.geometry().local(pair.first);
      auto& u = pair.second;
      const auto s = dynamic_cast<const EntropyFluxType*>(&analytical_flux_)
                         ->derived_local_function(entity)
                         ->get_alpha(x_in_inside_coords, u, param_)
                         .second;

      // if regularization was needed, we also need to replace u_n in that cell by its regularized version
      if (s > 0.) {
        std::cout << "regularization in limiter: time: " << param_.get("t")[0]
                  << ", entitycenter: " << XT::Common::to_string(entity.geometry().center())
                  << ", s: " << XT::Common::to_string(s, 15) << std::endl;
        const auto u_iso = dynamic_cast<const EntropyFluxType*>(&analytical_flux_)
                               ->basis_functions()
                               .calculate_isotropic_distribution(u)
                               .first;
        u = convex_combination(u, u_iso, s);
      } // if (s > 0)
    } // local_reconstructed_values
  } // void apply_limiter(...)

protected:
  RangeType convex_combination(const RangeType& u, const RangeType& u_bar, const RangeFieldType& theta)
  {
    RangeType u_scaled = u;
    u_scaled *= 1 - theta;
    RangeType u_bar_scaled = u_bar;
    u_bar_scaled *= theta;
    return u_scaled + u_bar_scaled;
  }

  const AnalyticalFluxType& analytical_flux_;
  const DiscreteFunctionType& source_;
  ReconstructedFunctionType& reconstructed_function_;
  const BasisfunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const RangeFieldType epsilon_;
  const std::vector<RangeType>& basis_values_;
  XT::Common::Parameter param_;
};

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp = int>
class NonLimitingLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::EntityType;

  template <class... Args>
  NonLimitingLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  void apply_local(const EntityType& /*entity*/)
  {
  }
}; // class NonLimitingLocalRealizabilityLimiter


template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class PositivityLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  template <class... Args>
  PositivityLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

    // vector to store thetas for each local reconstructed value
    thread_local std::vector<RangeFieldType> thetas(local_reconstructed_values.size());
    std::fill(thetas.begin(), thetas.end(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      for (size_t ii = 0; ii < u_l.size(); ++ii) {
        if (u_l[ii] >= u_bar[ii])
          continue;
        thetas[ll] = std::max(thetas[ll], u_l[ii] / (u_l[ii] - u_bar[ii]));
      }
    } // ll
    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    BaseType::apply_limiter(entity, theta_entity, local_reconstructed_values, u_bar);
  } // void apply_local(...)

private:
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::epsilon_;
}; // class PositivityLocalRealizabilityLimiter

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class DgLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::AnalyticalFluxType;
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

  template <class... Args>
  DgLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , triangulation_(basis_functions_.triangulation())
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
    BaseType::apply_limiter(entity, theta_entity, local_reconstructed_values, u_bar);
  } // void apply_local(...)

private:
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::basis_functions_;
  using BaseType::epsilon_;
  typename BasisfunctionImp::TriangulationType triangulation_;
}; // class DgLocalRealizabilityLimiter

#if HAVE_QHULL

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class ConvexHullLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  using typename BaseType::ReconstructedFunctionType;
  using typename BaseType::BasisfunctionType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  typedef typename std::vector<std::pair<RangeType, RangeFieldType>> PlaneCoefficientsType;

  template <class... Args>
  ConvexHullLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
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

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
bool ConvexHullLocalRealizabilityLimiter<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>::is_instantiated_ =
    false;

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
std::shared_ptr<typename ConvexHullLocalRealizabilityLimiter<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>::
                    PlaneCoefficientsType>
    ConvexHullLocalRealizabilityLimiter<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>::plane_coefficients_;

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class DgConvexHullLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using ThisType = DgConvexHullLocalRealizabilityLimiter;
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  using typename BaseType::ReconstructedFunctionType;
  using typename BaseType::BasisfunctionType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t block_size = (dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  typedef FieldVector<RangeFieldType, block_size - 1> BlockRangeType;
  typedef typename std::vector<std::pair<BlockRangeType, RangeFieldType>> BlockPlaneCoefficientsType;
  typedef FieldVector<BlockPlaneCoefficientsType, num_blocks> PlaneCoefficientsType;

public:
  template <class... Args>
  DgConvexHullLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
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
      const auto factor = basis_functions_.realizability_limiter_max(u_l, u_bar);
      u_l /= factor;
      auto u_bar_scaled = u_bar;
      u_bar_scaled /= factor;
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t offset = jj * block_size;
        if (u_l[offset] >= u_bar_scaled[offset])
          continue;
        thetas[ll] = std::max(thetas[ll], u_l[offset] / (u_l[offset] - u_bar_scaled[offset]));
      }
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t offset = jj * block_size;
        BlockRangeType q_k, q_bar_k;
        for (size_t mm = 0; mm < block_size - 1; ++mm) {
          q_k[mm] = u_l[offset + mm + 1] / u_l[offset];
          q_bar_k[mm] = u_bar_scaled[offset + mm + 1] / u_bar_scaled[offset];
        }
        for (const auto& coeffs : (*plane_coefficients_)[jj]) {
          const BlockRangeType& a = coeffs.first;
          const RangeFieldType& b = coeffs.second;
          RangeFieldType theta_li = (b - a * q_k) / (a * (q_bar_k - q_k));
          if (!(theta_li > 1.))
            thetas[ll] = std::max(thetas[ll], theta_li);
        } // coeffs
      } // jj
    } // ll

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    BaseType::apply_limiter(entity, theta_entity, local_reconstructed_values, u_bar);
  } // void apply_local(...)

private:
  // calculate half space representation of realizable set
  void calculate_plane_coefficients()
  {
    plane_coefficients_ = std::make_shared<PlaneCoefficientsType>();
    FieldVector<std::vector<FieldVector<RangeFieldType, block_size - 1>>, num_blocks> points;
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
        for (size_t ll = 0; ll < block_size - 1; ++ll)
          points[jj][ii][ll] = val[block_size * jj + 1 + ll];
      } // ii
      points[jj][blocked_quadrature[jj].size()] = FieldVector<RangeFieldType, block_size - 1>(0.);
    }
    std::vector<std::thread> threads(num_blocks);
    // Launch a group of threads
    for (size_t jj = 0; jj < num_blocks; ++jj)
      threads[jj] = std::thread(&ThisType::calculate_plane_coefficient_block, this, std::ref(points[jj]), jj);
    // Join the threads with the main thread
    for (size_t jj = 0; jj < num_blocks; ++jj)
      threads[jj].join();
  }

  void calculate_plane_coefficient_block(std::vector<FieldVector<RangeFieldType, block_size - 1>>& points, size_t jj)
  {
    orgQhull::Qhull qhull;
    qhull.runQhull("Realizable set", int(block_size) - 1, int(points.size()), &(points[0][0]), "Qt T1");
    const auto facet_end = qhull.endFacet();
    BlockPlaneCoefficientsType block_plane_coefficients(qhull.facetList().count());
    std::cout << "num_vertices: " << qhull.vertexList().count() << std::endl;
    size_t ii = 0;
    for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ii) {
      for (size_t ll = 0; ll < block_size - 1; ++ll)
        block_plane_coefficients[ii].first[ll] = *(facet.hyperplane().coordinates() + ll);
      block_plane_coefficients[ii].second = -facet.hyperplane().offset();
    } // ii
    (*plane_coefficients_)[jj] = block_plane_coefficients;
  }


  using BaseType::basis_functions_;
  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<PlaneCoefficientsType> plane_coefficients_;
}; // class ConvexHullLocalRealizabilityLimiter

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
bool DgConvexHullLocalRealizabilityLimiter<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>::is_instantiated_ =
    false;

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
std::shared_ptr<typename DgConvexHullLocalRealizabilityLimiter<AnalyticalFluxImp,
                                                               DiscreteFunctionImp,
                                                               BasisfunctionImp>::PlaneCoefficientsType>
    DgConvexHullLocalRealizabilityLimiter<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>::
        plane_coefficients_;

#else // HAVE_QHULL

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class ConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionImp>::value, "You are missing Qhull!");
};

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class DgConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionImp>::value, "You are missing Qhull!");
};

#endif // HAVE_QHULL

#if HAVE_CLP
template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class ClpLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  template <class... Args>
  ClpLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);
    assert(local_reconstructed_values.size() == 2 * dimDomain);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size());

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      // if (XT::Common::FloatCmp::eq(u_bar, u_l) || basis_functions_.obviously_realizable(u_l)) {
      if (XT::Common::FloatCmp::eq(u_bar, u_l)) {
        thetas[ll] = -epsilon_;
        continue;
      }

      // solve LP:
      // min \theta s.t.
      // (\sum x_i v_i) + \theta (u - \bar{u}) = u
      // x_i, \theta >= 0
      auto theta = solve_linear_program(u_bar, u_l);
      thetas[ll] = theta;
    } // ll

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    BaseType::apply_limiter(entity, theta_entity, local_reconstructed_values, u_bar);
  } // void apply_local(...)

private:
  void setup_linear_program()
  {
    if (!*lp_) {
      // We start with creating a model with dimRange rows and num_quad_points+1 columns */
      constexpr int num_rows = static_cast<int>(dimRange);
      assert(quadrature_.size() < std::numeric_limits<int>::max());
      int num_cols = static_cast<int>(quadrature_.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
      *lp_ = std::make_unique<ClpSimplex>(false);
      auto& lp = **lp_;
      // set number of rows
      lp.resize(num_rows, 0);

      // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
      // 0..num_rows
      std::array<int, num_rows> row_indices;
      for (size_t ii = 0; ii < num_rows; ++ii)
        row_indices[ii] = ii;

      // set columns for quadrature points
      assert(int(basis_values_.size()) == num_cols - 1);
      for (int ii = 0; ii < int(basis_values_.size()); ++ii) {
        const auto& v_i = basis_values_[ii];
        // First argument: number of elements in column
        // Second/Third argument: indices/values of column entries
        // Fourth/Fifth argument: lower/upper column bound, i.e. lower/upper bound for x_i. As all x_i should be
        // positive, set to 0/inf, which is the default.
        // Sixth argument: Prefactor in objective for x_i, this is 0 for all x_i, which is also the default;
        lp.addColumn(num_rows, row_indices.data(), &(v_i[0]));
      }

      // add theta column (set to random values, will be set correctly in solve_linear_program)
      // The bounds for theta should be [0,1], but we allow allow theta to be slightly
      // negative so we can check if u_l is on the boundary and if so, move it a little
      // away from the boundary
      // Also sets the prefactor in the objective to 1 for theta.
      lp.addColumn(num_rows, row_indices.data(), &(basis_values_[0][0]), -0.1, 1., 1.);
      lp.setLogLevel(0);
    } // if (!lp_)
  }

  //  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l, const size_t index)
  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l)
  {
    setup_linear_program();
    auto& lp = **lp_;
    constexpr int num_rows = static_cast<int>(dimRange);
    int num_cols = static_cast<int>(quadrature_.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
    RangeFieldType theta;
    const auto u_l_minus_u_bar = u_l - u_bar;

    // set rhs (equality constraints, so set both bounds equal
    for (size_t ii = 0; ii < num_rows; ++ii) {
      lp.setRowLower(ii, u_l[ii]);
      lp.setRowUpper(ii, u_l[ii]);
    }

    // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
    // 0..num_rows
    std::array<int, num_rows> row_indices;
    for (size_t ii = 0; ii < num_rows; ++ii)
      row_indices[ii] = ii;

    // delete and reset theta column
    const int last_col = num_cols - 1;
    lp.deleteColumns(1, &last_col);
    lp.addColumn(num_rows, row_indices.data(), &(u_l_minus_u_bar[0]), -0.1, 1., 1.);

    // Now solve
    lp.primal();
    theta = lp.objectiveValue();
    if (!lp.isProvenOptimal())
      theta = 1.;

    return theta;
  }

  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  using BaseType::basis_values_;
  XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>> lp_;
}; // class ClpLocalRealizabilityLimiter

#else // HAVE_CLP

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class ClpLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionImp>::value, "You are missing Clp!");
};

#endif // HAVE_CLP

template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class LPLocalRealizabilityLimiter
    : public LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>
{
  using BaseType = LocalRealizabilityLimiterBase<AnalyticalFluxImp, DiscreteFunctionImp, BasisfunctionImp>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::QuadratureType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;

  template <class... Args>
  LPLocalRealizabilityLimiter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);
    assert(local_reconstructed_values.size() == 2 * dimDomain);

    // get cell average
    const RangeType u_bar =
        source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size());

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      const auto& u_l = pair.second;
      // if (XT::Common::FloatCmp::eq(u_bar, u_l) || basis_functions_.obviously_realizable(u_l)) {
      if (XT::Common::FloatCmp::eq(u_bar, u_l)) {
        thetas[ll] = -epsilon_;
        continue;
      }

      // solve LP:
      // min \theta s.t.
      // (\sum x_i v_i) + \theta (u - \bar{u}) = u
      // x_i, \theta >= 0
      auto theta = solve_linear_program(u_bar, u_l);
      thetas[ll] = theta;
    } // ll

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    BaseType::apply_limiter(entity, theta_entity, local_reconstructed_values, u_bar);
  } // void apply_local(...)

private:
  void setup_linear_program()
  {
    if (!*lp_) {
      // We start with creating a model with dimRange rows and num_quad_points+1 columns */
      constexpr int num_rows = static_cast<int>(dimRange);
      assert(quadrature_.size() < std::numeric_limits<int>::max());
      int num_cols = static_cast<int>(quadrature_.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
      *lp_ = XT::Common::lp_solve::make_lp(num_rows, num_cols);
      auto& lp = **lp_;

      /* let us name our variables. Not required, but can be useful for debugging */
      for (int ii = 1; ii <= num_cols - 1; ++ii) {
        std::string ii_string = XT::Common::to_string(ii);
        std::string name_string = "x";
        for (size_t jj = 5 - ii_string.size(); jj > 0; --jj)
          name_string += "0";
        name_string += ii_string;
        XT::Common::lp_solve::set_col_name(lp, ii, name_string);
      }

      std::string name_string = "theta ";
      XT::Common::lp_solve::set_col_name(lp, num_cols, name_string);

      // In the call to set_column, the first entry (row 0) is the value of the objective function
      // (c_i in the objective function c^T x), the other entries are the entries of the i-th column
      // in the constraints matrix. The 0-th column is the rhs vector. The entry (0, 0) corresponds
      // to the initial value of the objective function.
      std::array<double, num_rows + 1> column;

      // set columns for quadrature points
      column[0] = 0.;
      assert(int(basis_values_.size()) == num_cols - 1);
      for (int ii = 0; ii < int(basis_values_.size()); ++ii) {
        const auto& v_i = basis_values_[ii];
        std::copy(v_i.begin(), v_i.end(), column.begin() + 1);
        XT::Common::lp_solve::set_column(lp, ii + 1, column.data());
      }

      // set all contraints to equality constraints
      for (int ii = 1; ii <= num_rows; ++ii)
        XT::Common::lp_solve::set_constr_type(lp, ii, XT::Common::lp_solve::eq());

      // Set bounds for variables. 0 <= x <= inf is the default for all variables.
      // The bounds for theta should be [0,1], but we allow allow theta to be slightly
      // negative so we can check if u_l is on the boundary and if so, move it a little
      // away from the boundary
      XT::Common::lp_solve::set_bounds(lp, num_cols, -0.1, 1.);

      /* I only want to see important messages on screen while solving */
      XT::Common::lp_solve::set_verbose(lp, XT::Common::lp_solve::important());
    } // if (!lp_)
  }

  //  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l, const size_t index)
  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l)
  {
    setup_linear_program();
    auto& lp = **lp_;
    constexpr int num_rows = static_cast<int>(dimRange);
    int num_cols = static_cast<int>(quadrature_.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
    RangeFieldType theta;
    const auto u_l_minus_u_bar = u_l - u_bar;

    // set rhs (column 0)
    std::array<double, num_rows + 1> column;
    column[0] = 0.;
    std::copy(u_l.begin(), u_l.end(), column.begin() + 1);
    XT::Common::lp_solve::set_rh_vec(lp, column.data());
    XT::Common::lp_solve::set_rh(lp, 0, column[0]);

    // set theta column
    column[0] = 1.;
    std::copy(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end(), column.begin() + 1);
    XT::Common::lp_solve::set_column(lp, num_cols, column.data());

    /* Now let lpsolve calculate a solution */
    const auto solve_status = XT::Common::lp_solve::solve(lp);
    theta = XT::Common::lp_solve::get_objective(lp);
    if (solve_status != XT::Common::lp_solve::optimal()) {
      //      XT::Common::lp_solve::write_LP(lp, stdout);
      theta = 1.;
    }

    return theta;
  }

  using BaseType::source_;
  using BaseType::reconstructed_function_;
  using BaseType::quadrature_;
  using BaseType::epsilon_;
  using BaseType::basis_values_;
  XT::Common::PerThreadValue<std::unique_ptr<typename XT::Common::lp_solve::LinearProgram>> lp_;
}; // class LPLocalRealizabilityLimiter


template <class LocalRealizabilityLimiterImp, class Traits>
class RealizabilityLimiter;


namespace internal {


template <class LocalRealizabilityLimiterImp>
struct RealizabilityLimiterTraits
{
  using LocalRealizabilityLimiterType = LocalRealizabilityLimiterImp;
  using AnalyticalFluxType = typename LocalRealizabilityLimiterImp::AnalyticalFluxType;
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
  using AnalyticalFluxType = typename Traits::AnalyticalFluxType;
  using BasisfunctionType = typename Traits::BasisfunctionType;
  using QuadratureType = typename Traits::QuadratureType;
  using RangeFieldType = typename Traits::RangeFieldType;
  using RangeType = typename Traits::RangeType;
  using ReconstructedFunctionType = typename Traits::ReconstructedFunctionType;

  RealizabilityLimiter(const AnalyticalFluxType& analytical_flux,
                       const BasisfunctionType& basis_functions,
                       const QuadratureType& quadrature,
                       const RangeFieldType epsilon = 1e-10)
    : analytical_flux_(analytical_flux)
    , basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , epsilon_(epsilon)
    , basis_values_(quadrature_.size())
  {
    for (size_t ii = 0; ii < quadrature_.size(); ++ii)
      basis_values_[ii] = basis_functions_.evaluate(quadrature_[ii].position());
  }

  template <class SourceType>
  void apply(const SourceType& source, ReconstructedFunctionType& range, const XT::Common::Parameter& param) const
  {
    static_assert(is_discrete_function<SourceType>::value,
                  "SourceType has to be derived from DiscreteFunction (use the non-reconstructed values!)");
    LocalRealizabilityLimiterType local_realizability_limiter(
        analytical_flux_, source, range, basis_functions_, quadrature_, epsilon_, basis_values_, param);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(source.space().grid_layer());
    walker.append(local_realizability_limiter);
    walker.walk(true);
  } // void apply(...)

private:
  const AnalyticalFluxType& analytical_flux_;
  const BasisfunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const RangeFieldType epsilon_;
  std::vector<RangeType> basis_values_;
}; // class RealizabilityLimiter<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
