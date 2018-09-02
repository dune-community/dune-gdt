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

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>

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

#include <dune/gdt/operators/fv/reconstruction/reconstructed_function.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/local/fluxes/entropybased.hh>

namespace Dune {
namespace GDT {


// The limiters in this file act on the reconstructed values on the reconstructed values to ensure realizability.
// Another possibility is to include the realizability limiting in the slope limiting, see
// dune/gdt/operators/fv/reconstruction/slopes.hh.
// Which one of the two limiting approaches is more suitable depends on the actual reconstruction. If the values are
// reconstructed on the midpoint of each intersection only, the slope limiting is easy to use and it is possible to
// limit in characteristic variables at no extra computational cost. Moreover, if the the values are not only
// reconstructed at the interfaces but also within the cell, limiting the slope before calculating the reconstructed
// values might be less expensive than limiting the reconstructed values afterwards. However, if we have a quadrature
// with several points on each interface, we do a dimension by dimension reconstruction on the whole stencil, and in the
// current implementation slope limiting is performed at every step. It might be enough to perform the realizability
// limiting in the last step only, which would bring the slope limiters on a par with the reconstructions in this file
// performance-wise.
template <class AnalyticalFluxImp, class DiscreteFunctionImp, class BasisfunctionImp>
class LocalRealizabilityLimiterBase
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionImp::SpaceType::GridLayerType>
{
public:
  using AnalyticalFluxType = AnalyticalFluxImp;
  using DiscreteFunctionType = DiscreteFunctionImp;
  using BasisfunctionType = BasisfunctionImp;
  using GridLayerType = typename DiscreteFunctionType::SpaceType::GridLayerType;
  using EntityType = typename DiscreteFunctionType::EntityType;
  using DomainType = typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate;
  using RangeType = typename DiscreteFunctionType::RangeType;
  using IndexSetType = typename GridLayerType::IndexSet;
  using RangeFieldType = typename DiscreteFunctionType::RangeFieldType;
  using DomainFieldType = typename DiscreteFunctionType::DomainFieldType;
  static const size_t dimDomain = DiscreteFunctionType::dimDomain;
  static const size_t dimRange = DiscreteFunctionType::dimRange;
  using ReconstructedFunctionType =
      ReconstructedLocalizableFunction<GridLayerType, DomainFieldType, dimDomain, RangeFieldType, dimRange>;
  using EntropyFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, DiscreteFunctionType>;
  using LocalReconstructedValuesType = typename ReconstructedFunctionType::LocalReconstructedValuesType;

  LocalRealizabilityLimiterBase(const AnalyticalFluxType& analytical_flux,
                                const DiscreteFunctionType& source,
                                ReconstructedFunctionType& reconstructed_function,
                                const BasisfunctionType& basis_functions,
                                const RangeFieldType epsilon,
                                const std::vector<RangeType>& basis_values,
                                const XT::Common::Parameter& param)
    : analytical_flux_(analytical_flux)
    , source_(source)
    , reconstructed_function_(reconstructed_function)
    , basis_functions_(basis_functions)
    , epsilon_(epsilon)
    , basis_values_(basis_values)
    , param_(param)
  {
    param_.set("boundary", {0.});
  }

  virtual ~LocalRealizabilityLimiterBase()
  {
  }

  // scalar version
  void apply_limiter(const EntityType& entity,
                     const RangeFieldType theta_entity,
                     LocalReconstructedValuesType& local_reconstructed_values,
                     const RangeType& u_bar,
                     bool add_epsilon = true)
  {
    // apply slope limiter
    assert(dynamic_cast<const EntropyFluxType*>(&analytical_flux_) != nullptr
           && "analytical_flux_ has to be derived from EntropyBasedLocalFlux");
    auto theta = add_epsilon ? theta_entity + epsilon_ : theta_entity;
    if (theta > 0.) {
      // std::cout << "limited with theta: " << theta << " and epsilon " << epsilon_ << std::endl;
      if (theta > 1.)
        theta = 1.;
      for (auto& pair : local_reconstructed_values) {
        auto& u = pair.second;
        u = convex_combination(u, u_bar, theta);
      }
    }
    check_solvability(entity, local_reconstructed_values, u_bar);
  } // void apply_limiter(...)

  // version for component-wise limiting
  void apply_limiter(const EntityType& entity,
                     const RangeType& theta_entity,
                     LocalReconstructedValuesType& local_reconstructed_values,
                     const RangeType& u_bar,
                     bool add_epsilon = true)
  {
    // apply slope limiter
    assert(dynamic_cast<const EntropyFluxType*>(&analytical_flux_) != nullptr
           && "analytical_flux_ has to be derived from EntropyBasedLocalFlux");
    for (size_t ii = 0; ii < dimRange; ++ii) {
      auto theta_ii = add_epsilon ? theta_entity[ii] + epsilon_ : theta_entity[ii];
      if (theta_ii > 0.) {
        // std::cout << "limited with theta: " << theta_ii << " and epsilon " << epsilon_ << std::endl;
        if (theta_ii > 1.)
          theta_ii = 1.;
        for (auto& pair : local_reconstructed_values) {
          auto& u_ii = pair.second[ii];
          u_ii = convex_combination(u_ii, u_bar[ii], theta_ii);
        }
      } // if (theta_ii > 0)
    } // ii
    check_solvability(entity, local_reconstructed_values, u_bar);
  }

protected:
  void check_solvability(const EntityType& entity,
                         LocalReconstructedValuesType& local_reconstructed_values,
                         const RangeType& u_bar)

  {
    // Try to solve optimization problem for all reconstructed values. If it fails for one value,
    // we cannot guarantee realizability preservation, so disable reconstruction in that case.
    const auto local_func = dynamic_cast<const EntropyFluxType*>(&analytical_flux_)->derived_local_function(entity);
    std::pair<DomainType, RangeType> current_pair = *local_reconstructed_values.begin();
    try {
      for (const auto& pair : local_reconstructed_values) {
        current_pair = pair;
        const auto x_in_inside_coords = entity.geometry().local(pair.first);
        const auto& u = pair.second;
        local_func->get_alpha(x_in_inside_coords, u, param_, false, false);
      } // local_reconstructed_values
    } catch (const Dune::MathError&) {
      std::cout << "Reconstruction disabled at time " << XT::Common::to_string(param_.get("t")[0], 15)
                << " and entity with center " << XT::Common::to_string(entity.geometry().center(), 15) << std::endl;
      std::cout << "Solving failed for moments " << XT::Common::to_string(current_pair.second, 15)
                << " at x = " << XT::Common::to_string(current_pair.first, 15) << std::endl;
      // solving failed for reconstructed value, so check that it works with u_bar ...
      local_func->get_alpha(entity.geometry().local(entity.geometry().center()), u_bar, param_, false, false);
      // ... and set all reconstructed values to u_bar
      for (auto& pair : local_reconstructed_values)
        pair.second = u_bar;
    }
  }

  RangeType convex_combination(const RangeType& u, const RangeType& u_bar, const RangeFieldType& theta)
  {
    RangeType u_scaled = u;
    u_scaled *= 1 - theta;
    RangeType u_bar_scaled = u_bar;
    u_bar_scaled *= theta;
    return u_scaled + u_bar_scaled;
  }

  RangeFieldType convex_combination(const RangeFieldType& u, const RangeFieldType& u_bar, const RangeFieldType& theta)
  {
    return u_bar * theta + u * (1. - theta);
  }

  const AnalyticalFluxType& analytical_flux_;
  const DiscreteFunctionType& source_;
  ReconstructedFunctionType& reconstructed_function_;
  const BasisfunctionType& basis_functions_;
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

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);
    const auto u_bar = source_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));
    BaseType::apply_limiter(entity, 0., local_reconstructed_values, u_bar, false);
  }

private:
  using BaseType::reconstructed_function_;
  using BaseType::source_;
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
    RangeType thetas(0.);

    for (const auto& pair : local_reconstructed_values) {
      const auto& u = pair.second;
      for (size_t ii = 0; ii < u.size(); ++ii) {
        if (u[ii] >= u_bar[ii])
          continue;
        if (u_bar[ii] < epsilon_)
          thetas[ii] = 1.;
        else if (u[ii] < epsilon_)
          thetas[ii] = std::max(thetas[ii], (epsilon_ - u[ii]) / (u_bar[ii] - u[ii]));
      } // ii
    } // ll

    const auto theta_max = *std::max_element(thetas.begin(), thetas.end());
    BaseType::apply_limiter(entity, theta_max, local_reconstructed_values, u_bar, false);
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
    RangeType thetas(0.);

    for (const auto& pair : local_reconstructed_values) {
      const auto& u = pair.second;
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        const auto& u0 = u[2 * ii];
        const auto& u1 = u[2 * ii + 1];
        const auto& ubar0 = u_bar[2 * ii];
        const auto& ubar1 = u_bar[2 * ii + 1];
        const auto& vj = triangulation_[ii];
        const auto& vjplus1 = triangulation_[ii + 1];
        FieldVector<RangeFieldType, 3> thetas_ii;
        if (!is_epsilon_realizable(ubar0, ubar1, vj, vjplus1, epsilon_)) {
          thetas[2 * ii] = 1.;
        } else {
          thetas_ii[0] = (epsilon_ - u0) / (ubar0 - u0);
          thetas_ii[1] =
              (u0 * vj - u1 + epsilon_ * std::sqrt(std::pow(vj, 2) + 1)) / ((ubar1 - u1) - (ubar0 - u0) * vj);
          thetas_ii[2] = (u0 * vjplus1 - u1 - epsilon_ * std::sqrt(std::pow(vjplus1, 2) + 1))
                         / ((ubar1 - u1) - (ubar0 - u0) * vjplus1);
          for (size_t kk = 0; kk < 3; ++kk)
            if (thetas_ii[kk] >= 0. && thetas_ii[kk] <= 1.)
              thetas[2 * ii] = std::max(thetas[2 * ii], thetas_ii[kk]);
        } // else (!realizable)
        thetas[2 * ii + 1] = thetas[2 * ii];
      } // ii
    } // local_reconstructed_values
    BaseType::apply_limiter(entity, thetas, local_reconstructed_values, u_bar, false);
  } // void apply_local(...)

private:
  bool is_epsilon_realizable(const RangeFieldType ubar0,
                             const RangeFieldType ubar1,
                             const RangeFieldType v0,
                             const RangeFieldType v1,
                             const RangeFieldType eps) const
  {
    bool ret = (ubar0 >= eps) && (ubar1 <= v1 * ubar0 - eps * std::sqrt(std::pow(v1, 2) + 1))
               && (v0 * ubar0 + eps * std::sqrt(std::pow(v0, 2) + 1) <= ubar1);
    return ret;
  }

  using BaseType::basis_functions_;
  using BaseType::epsilon_;
  using BaseType::reconstructed_function_;
  using BaseType::source_;
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

    size_t ll = static_cast<size_t>(-1);
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      // rescale u_l, u_bar
      auto u_l = pair.second;
      auto u_bar_minus_u_l = u_bar - u_l;
      const auto factor = std::max(basis_functions_.density(u_l), basis_functions_.density(u_bar)) / 2.;
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
    const auto& quadrature = basis_functions_.quadratures().merged();
    std::vector<FieldVector<RangeFieldType, dimRange>> points(quadrature.size() + 1);
    points[0] = FieldVector<RangeFieldType, dimRange>(0);
    size_t ii = 1;
    for (const auto& quad_point : quadrature)
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
  using typename BaseType::ReconstructedFunctionType;
  using typename BaseType::BasisfunctionType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t block_size = (dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  typedef FieldVector<RangeFieldType, block_size> BlockRangeType;
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

    size_t ll = static_cast<size_t>(-1);
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      // rescale u_l, u_bar
      auto u_l = pair.second;
      if (XT::Common::FloatCmp::eq(u_l, u_bar))
        continue;
      const auto factor = std::max(basis_functions_.density(u_l), basis_functions_.density(u_bar)) / 2.;
      u_l /= factor;
      auto u_bar_scaled = u_bar;
      u_bar_scaled /= factor;
      // check positivity of first moment on each spherical triangle
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t offset = jj * block_size;
        if (u_l[offset] >= u_bar_scaled[offset])
          continue;
        thetas[ll] = std::max(thetas[ll], u_l[offset] / (u_l[offset] - u_bar_scaled[offset]));
      }
      // check convex hull property of other moments
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t offset = jj * block_size;
        BlockRangeType q_k, q_bar_k;
        for (size_t mm = 0; mm < block_size - 1; ++mm) {
          q_k[mm] = u_l[offset + mm + 1];
          q_bar_k[mm] = u_bar_scaled[offset + mm + 1];
        }
        for (const auto& coeffs : (*plane_coefficients_)[jj]) {
          const BlockRangeType& a = coeffs.first;
          const RangeFieldType& b = coeffs.second;
          const auto u0 = u_l[offset];
          const auto ubar0 = u_bar_scaled[offset];
          RangeFieldType theta_li = (b * u0 - a * q_k) / (a * (q_bar_k - q_k) - b * (ubar0 - u0));
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
    const auto& quadratures = basis_functions_.quadratures();
    plane_coefficients_ = std::make_shared<PlaneCoefficientsType>();
    FieldVector<std::vector<FieldVector<RangeFieldType, block_size>>, num_blocks> points;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      points[jj].resize(quadratures[jj].size() + 1);
      for (size_t ii = 0; ii < quadratures[jj].size(); ++ii) {
        const auto val = basis_functions_.evaluate(quadratures()[jj][ii].position(), jj);
        for (size_t ll = 0; ll < block_size; ++ll)
          points[jj][ii][ll] = val[block_size * jj + ll];
      } // ii
      points[jj][quadratures[jj].size()] = FieldVector<RangeFieldType, block_size>(0.);
    }
    std::vector<std::thread> threads(num_blocks);
    // Launch a group of threads
    for (size_t jj = 0; jj < num_blocks; ++jj)
      threads[jj] = std::thread(&ThisType::calculate_plane_coefficient_block, this, std::ref(points[jj]), jj);
    // Join the threads with the main thread
    for (size_t jj = 0; jj < num_blocks; ++jj)
      threads[jj].join();
  }

  void calculate_plane_coefficient_block(std::vector<FieldVector<RangeFieldType, block_size>>& points, size_t jj)
  {
    orgQhull::Qhull qhull;
    boost::iostreams::stream<boost::iostreams::null_sink> null_ostream((boost::iostreams::null_sink()));
    qhull.setOutputStream(&null_ostream);
    qhull.setErrorStream(&null_ostream);
    qhull.runQhull("Realizable set", int(block_size), int(points.size()), &(points[0][0]), "Qt T1");
    const auto facet_end = qhull.endFacet();
    BlockPlaneCoefficientsType block_plane_coefficients(qhull.facetList().count());
    //    std::cout << "num_vertices: " << qhull.vertexList().count() << std::endl;
    size_t ii = 0;
    for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ii) {
      for (size_t ll = 0; ll < block_size; ++ll)
        block_plane_coefficients[ii].first[ll] = *(facet.hyperplane().coordinates() + ll);
      block_plane_coefficients[ii].second = -facet.hyperplane().offset();
    } // ii
    (*plane_coefficients_)[jj] = block_plane_coefficients;
  }

  using BaseType::basis_functions_;
  using BaseType::source_;
  using BaseType::reconstructed_function_;
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

    RangeFieldType theta_entity(-epsilon_);

    for (const auto& pair : local_reconstructed_values) {
      const auto& u = pair.second;
      if (XT::Common::FloatCmp::eq(u_bar, u))
        continue;

      // solve LP:
      // min \theta s.t.
      // (\sum x_i v_i) + \theta (u - \bar{u}) = u
      // x_i, \theta >= 0
      theta_entity = std::max(theta_entity, solve_linear_program(u_bar, u));
    } // ll
    BaseType::apply_limiter(entity, theta_entity, local_reconstructed_values, u_bar, true);
  } // void apply_local(...)

private:
  void setup_linear_program()
  {
    if (!*lp_) {
      const auto& quadrature = basis_functions_.quadratures().merged();
      // We start with creating a model with dimRange rows and num_quad_points+1 columns */
      constexpr int num_rows = static_cast<int>(dimRange);
      assert(quadrature.size() < std::numeric_limits<int>::max());
      int num_cols = static_cast<int>(quadrature.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
      *lp_ = std::make_unique<ClpSimplex>(false);
      auto& lp = **lp_;
      // set number of rows
      lp.resize(num_rows, 0);

      // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
      // 0..num_rows
      std::array<int, num_rows> row_indices;
      for (int ii = 0; ii < num_rows; ++ii)
        row_indices[ii] = ii;

      // set columns for quadrature points
      assert(int(basis_values_.size()) == num_cols - 1);
      for (int ii = 0; ii < num_cols - 1; ++ii) {
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
    int num_cols = static_cast<int>(basis_functions_.quadratures().merged().size()
                                    + 1); /* variables are x_1, ..., x_{num_quad_points}, theta */
    RangeFieldType theta;
    const auto u_l_minus_u_bar = u_l - u_bar;

    // set rhs (equality constraints, so set both bounds equal
    for (int ii = 0; ii < num_rows; ++ii) {
      lp.setRowLower(ii, u_l[ii]);
      lp.setRowUpper(ii, u_l[ii]);
    }

    // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
    // 0..num_rows
    std::array<int, num_rows> row_indices;
    for (int ii = 0; ii < num_rows; ++ii)
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
  using BaseType::basis_functions_;
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


template <class LocalRealizabilityLimiterImp, class Traits>
class RealizabilityLimiter;


namespace internal {


template <class LocalRealizabilityLimiterImp>
struct RealizabilityLimiterTraits
{
  using LocalRealizabilityLimiterType = LocalRealizabilityLimiterImp;
  using AnalyticalFluxType = typename LocalRealizabilityLimiterImp::AnalyticalFluxType;
  using BasisfunctionType = typename LocalRealizabilityLimiterImp::BasisfunctionType;
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
  using RangeFieldType = typename Traits::RangeFieldType;
  using RangeType = typename Traits::RangeType;
  using ReconstructedFunctionType = typename Traits::ReconstructedFunctionType;

  RealizabilityLimiter(const AnalyticalFluxType& analytical_flux,
                       const BasisfunctionType& basis_functions,
                       const RangeFieldType epsilon = 1e-8)
    : analytical_flux_(analytical_flux)
    , basis_functions_(basis_functions)
    , epsilon_(epsilon)
    , basis_values_(basis_functions_.quadratures().merged().size())
  {
    const auto& quadrature = basis_functions_.quadratures().merged();
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      basis_values_[ii] = basis_functions_.evaluate(quadrature[ii].position());
  }

  template <class SourceType>
  void apply(const SourceType& source, ReconstructedFunctionType& range, const XT::Common::Parameter& param) const
  {
    static_assert(is_discrete_function<SourceType>::value,
                  "SourceType has to be derived from DiscreteFunction (use the non-reconstructed values!)");
    LocalRealizabilityLimiterType local_realizability_limiter(
        analytical_flux_, source, range, basis_functions_, epsilon_, basis_values_, param);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(source.space().grid_layer());
    walker.append(local_realizability_limiter);
    walker.walk(true);
  } // void apply(...)

private:
  const AnalyticalFluxType& analytical_flux_;
  const BasisfunctionType& basis_functions_;
  const RangeFieldType epsilon_;
  std::vector<RangeType> basis_values_;
}; // class RealizabilityLimiter<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
