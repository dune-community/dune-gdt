// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_SLOPES_HH
#define DUNE_GDT_OPERATORS_FV_SLOPES_HH

#if HAVE_CLP
#  include <coin/ClpSimplex.hpp>
#endif // HAVE_CLP

#if HAVE_QHULL
#  include <dune/xt/common/disable_warnings.hh>
#  include <libqhullcpp/Qhull.h>
#  include <libqhullcpp/QhullFacetList.h>
#  include <dune/xt/common/reenable_warnings.hh>
#endif // HAVE_QHULL

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parallel/threadstorage.hh>

#include <dune/gdt/momentmodels/basisfunctions/partial_moments.hh>
#include <dune/gdt/momentmodels/entropyflux.hh>

#include "internal.hh"

namespace Dune {
namespace GDT {


template <class E, class EigenVectorWrapperType, size_t stencil_size = 3>
class SlopeBase
{
public:
  virtual ~SlopeBase() = default;

  using VectorType = typename EigenVectorWrapperType::VectorType;
  using MatrixType = typename EigenVectorWrapperType::MatrixType;
  using StencilType = FieldVector<VectorType, stencil_size>;

  virtual SlopeBase* copy() const = 0;

  // at least one of the following two methods has to be implemented
  // returns (limited) slope in ordinary coordinates
  virtual VectorType
  get(const E& entity, const StencilType& stencil, const EigenVectorWrapperType& eigenvectors, const size_t dd) const
  {
    VectorType slope_char = get_char(entity, stencil, eigenvectors, dd);
    VectorType slope;
    eigenvectors.apply_eigenvectors(dd, slope_char, slope);
    return slope;
  }

  // returns (limited) slope in characteristic coordinates
  virtual VectorType get_char(const E& entity,
                              const StencilType& stencil,
                              const EigenVectorWrapperType& eigenvectors,
                              const size_t dd) const
  {
    VectorType slope = get(entity, stencil, eigenvectors, dd);
    // convert to characteristic coordinates
    VectorType slope_char;
    eigenvectors.apply_inverse_eigenvectors(dd, slope, slope_char);
    return slope_char;
  }
};


// Central slope without any limiting
template <class E, class EigenVectorWrapperType, size_t stencil_size>
class CentralSlope : public SlopeBase<E, EigenVectorWrapperType, stencil_size>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, stencil_size>;
  using ThisType = CentralSlope;

public:
  using typename BaseType::StencilType;
  using typename BaseType::VectorType;

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& /*entity*/,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& /*eigenvectors*/,
                         const size_t /*dd*/) const override final
  {
    const VectorType& u_left = stencil[0];
    const VectorType& u_right = stencil[2];
    return (u_right - u_left) / 2.;
  }
};


// Left slope without any limiting
template <class E, class EigenVectorWrapperType, size_t stencil_size>
class LeftSlope : public SlopeBase<E, EigenVectorWrapperType, stencil_size>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, stencil_size>;
  using ThisType = LeftSlope;

public:
  using typename BaseType::StencilType;
  using typename BaseType::VectorType;

  virtual BaseType* copy() const override final
  {
    return new ThisType();
  }

  virtual VectorType get(const E& /*entity*/,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& /*eigenvectors*/,
                         const size_t /*dd*/) const override final
  {
    const VectorType& u_left = stencil[0];
    const VectorType& u = stencil[1];
    return u - u_left;
  }
};


// Right slope without any limiting
template <class E, class EigenVectorWrapperType, size_t stencil_size>
class RightSlope : public SlopeBase<E, EigenVectorWrapperType, stencil_size>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, stencil_size>;
  using ThisType = RightSlope;

public:
  using typename BaseType::StencilType;
  using typename BaseType::VectorType;

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& /*entity*/,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& /*eigenvectors*/,
                         const size_t /*dd*/) const override final
  {
    const VectorType& u = stencil[1];
    const VectorType& u_right = stencil[2];
    return u_right - u;
  }
};


// Zero slope
template <class E, class EigenVectorWrapperType, size_t stencil_size>
class NoSlope : public SlopeBase<E, EigenVectorWrapperType, stencil_size>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, stencil_size>;
  using ThisType = NoSlope;
  using typename BaseType::VectorType;

public:
  using typename BaseType::StencilType;

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& /*entity*/,
                         const StencilType& /*stencil*/,
                         const EigenVectorWrapperType& /*eigenvectors*/,
                         const size_t /*dd*/) const override final
  {
    return VectorType(0.);
  }
};


template <class E, class EigenVectorWrapperType>
class MinmodSlope : public SlopeBase<E, EigenVectorWrapperType, 3>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using ThisType = MinmodSlope;
  using typename BaseType::VectorType;

public:
  using typename BaseType::StencilType;

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get_char(const E& /*entity*/,
                              const StencilType& stencil,
                              const EigenVectorWrapperType& eigenvectors,
                              const size_t dd) const override final
  {
    const VectorType slope_left = stencil[1] - stencil[0];
    const VectorType slope_right = stencil[2] - stencil[1];
    VectorType slope_left_char, slope_right_char;
    eigenvectors.apply_inverse_eigenvectors(dd, slope_left, slope_left_char);
    eigenvectors.apply_inverse_eigenvectors(dd, slope_right, slope_right_char);
    return minmod(slope_left_char, slope_right_char);
  }

  static VectorType minmod(const VectorType& first_slope, const VectorType& second_slope)
  {
    VectorType ret(0.);
    for (size_t ii = 0; ii < first_slope.size(); ++ii)
      if (std::signbit(first_slope[ii]) == std::signbit(second_slope[ii])) // check for equal sign
        ret[ii] = (std::abs(first_slope[ii]) < std::abs(second_slope[ii])) ? first_slope[ii] : second_slope[ii];
    return ret;
  }
};


template <class E, class EigenVectorWrapperType>
class McSlope : public SlopeBase<E, EigenVectorWrapperType, 3>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using ThisType = McSlope;
  using MinmodType = MinmodSlope<E, EigenVectorWrapperType>;

public:
  using typename BaseType::StencilType;
  using typename BaseType::VectorType;

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get_char(const E& /*entity*/,
                              const StencilType& stencil,
                              const EigenVectorWrapperType& eigenvectors,
                              const size_t dd) const override final
  {
    const VectorType& u_left = stencil[0];
    const VectorType& u = stencil[1];
    const VectorType& u_right = stencil[2];
    const VectorType slope_left = (u - u_left) * 2.;
    const VectorType slope_right = (u_right - u) * 2.;
    const VectorType slope_center = (u_right - u_left) / 2.;
    VectorType slope_left_char, slope_right_char, slope_center_char;
    eigenvectors.apply_inverse_eigenvectors(dd, slope_left, slope_left_char);
    eigenvectors.apply_inverse_eigenvectors(dd, slope_right, slope_right_char);
    eigenvectors.apply_inverse_eigenvectors(dd, slope_center, slope_center_char);
    const VectorType first_slope = MinmodType::minmod(slope_left_char, slope_right_char);
    return MinmodType::minmod(first_slope, slope_center_char);
  }
};


template <class E, class EigenVectorWrapperType>
class SuperbeeSlope : public SlopeBase<E, EigenVectorWrapperType, 3>
{
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using ThisType = SuperbeeSlope;
  using MinmodType = MinmodSlope<E, EigenVectorWrapperType>;

public:
  using typename BaseType::StencilType;
  using typename BaseType::VectorType;

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get_char(const E& /*entity*/,
                              const StencilType& stencil,
                              const EigenVectorWrapperType& eigenvectors,
                              const size_t dd) const override final
  {
    const VectorType slope_left = stencil[1] - stencil[0];
    const VectorType slope_right = stencil[2] - stencil[1];
    VectorType slope_left_char, slope_right_char;
    eigenvectors.apply_inverse_eigenvectors(dd, slope_left, slope_left_char);
    eigenvectors.apply_inverse_eigenvectors(dd, slope_right, slope_right_char);
    const VectorType first_slope = MinmodType::minmod(slope_left_char, slope_right_char * 2.);
    const VectorType second_slope = MinmodType::minmod(slope_left_char * 2., slope_right_char);
    return maxmod(first_slope, second_slope);
  }

  static VectorType maxmod(const VectorType& first_slope, const VectorType& second_slope)
  {
    VectorType ret(0.);
    for (size_t ii = 0; ii < first_slope.size(); ++ii)
      if (std::signbit(first_slope[ii]) == std::signbit(second_slope[ii])) // check for equal sign
        ret[ii] = (std::abs(first_slope[ii]) < std::abs(second_slope[ii])) ? second_slope[ii] : first_slope[ii];
    return ret;
  }
}; // class SuperbeeSlope


template <class GV, class MomentBasis>
class RealizabilityLimiterBase
{
public:
  using E = XT::Grid::extract_entity_t<GV>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  using StateType = typename EntropyFluxType::StateType;

  RealizabilityLimiterBase(const EntropyFluxType& entropy_flux)
    : entropy_flux_(entropy_flux)
    , local_flux_(entropy_flux_.derived_local_function())
  {}

  RealizabilityLimiterBase(const RealizabilityLimiterBase& other)
    : entropy_flux_(other.entropy_flux_)
    , local_flux_(entropy_flux_.derived_local_function())
  {}

  // Ensures we are able to solve the optimization problems for the reconstructed values without regularization. If
  // optimization fails, we just set the slope to 0.
  template <class VectorType>
  void ensure_solvability(const E& entity, const VectorType& u, VectorType& slope) const
  {
    local_flux_->bind(entity);
    try {
      const auto u_left = u - slope * 0.5;
      const auto u_right = u + slope * 0.5;
      local_flux_->get_alpha(u_left, false);
      local_flux_->get_alpha(u_right, false);
    } catch (const Dune::MathError&) {
      std::fill(slope.begin(), slope.end(), 0.);
    }
  }

private:
  const EntropyFluxType& entropy_flux_;
  std::unique_ptr<typename EntropyFluxType::Localfunction> local_flux_;
};


// Realizability limiter that ensures positivity of the components of u in noncharacteristic variables. Uses single
// limiter variable for all components.
template <class GV,
          class MomentBasis,
          class EigenVectorWrapperType,
          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>>
class PositivityLimitedSlope
  : public SlopeBase<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType, 3>
  , public RealizabilityLimiterBase<GV, MomentBasis>
{
  using ThisType = PositivityLimitedSlope;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  static const size_t dimRange = MomentBasis::dimRange;
  using RealizabilityBaseType = RealizabilityLimiterBase<GV, MomentBasis>;
  using typename RealizabilityBaseType::E;
  using typename RealizabilityBaseType::EntropyFluxType;
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using typename BaseType::VectorType;

public:
  using typename BaseType::StencilType;

  // This limiter ensures u_i >= epsilon for all components u_i of u.
  PositivityLimitedSlope(const EntropyFluxType& entropy_flux, const RangeFieldType epsilon = 0.)
    : RealizabilityBaseType(entropy_flux)
    , epsilon_(epsilon)
  {}

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& entity,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& eigenvectors,
                         const size_t dd) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType& u_bar = stencil[1];
    if (!is_epsilon_realizable(u_bar, epsilon_))
      return zero_vector;
    const VectorType slope = slope_limiter_.get(entity, stencil, eigenvectors, dd);
    // this needs to be changed for other interface quadratures (see above)
    const FieldVector<VectorType, 2> reconstructed_values{u_bar - 0.5 * slope, u_bar + 0.5 * slope};
    VectorType thetas(0.);
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u = reconstructed_values[kk];
      for (size_t ii = 0; ii < u.size(); ++ii) {
        if (u[ii] >= u_bar[ii])
          continue;
        else if (u[ii] < epsilon_)
          thetas[ii] = std::max(thetas[ii], (epsilon_ - u[ii]) / (u_bar[ii] - u[ii]));
      }
    } // kk

    const auto theta_max = *std::max_element(thetas.begin(), thetas.end());
    assert(XT::Common::FloatCmp::le(theta_max, 1.) && XT::Common::FloatCmp::ge(theta_max, 0.));
    VectorType ret = slope * (1 - theta_max);
    this->ensure_solvability(entity, u_bar, ret);
    return ret;
  }

private:
  static bool is_epsilon_realizable(const VectorType& u, const RangeFieldType epsilon)
  {
    for (const auto& val : u)
      if (val < epsilon)
        return false;
    return true;
  }

  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
}; // class PositivityLimitedSlope<...>


// Realizability limiter that ensures positivity of the components of u in noncharacteristic variables. Uses single
// limiter variable for all components.
template <class GV,
          class RangeFieldType,
          size_t dimRange,
          class EigenVectorWrapperType,
          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>>
class Dg1dRealizabilityLimitedSlope
  : public SlopeBase<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType, 3>
  , public RealizabilityLimiterBase<GV, PartialMomentBasis<typename GV::ctype, 1, RangeFieldType, dimRange, 1, 1>>
{
  using ThisType = Dg1dRealizabilityLimitedSlope;
  using MomentBasis = Dune::GDT::PartialMomentBasis<RangeFieldType, 1, RangeFieldType, dimRange, 1, 1, 1>;
  using RealizabilityBaseType = RealizabilityLimiterBase<GV, MomentBasis>;
  using typename RealizabilityBaseType::E;
  using typename RealizabilityBaseType::EntropyFluxType;
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using typename BaseType::VectorType;

public:
  using typename BaseType::StencilType;

  // This limiter ensures u_i >= epsilon for all components u_i of u.
  Dg1dRealizabilityLimitedSlope(const EntropyFluxType& entropy_flux,
                                const MomentBasis& basis_functions,
                                const RangeFieldType epsilon = 0.)
    : RealizabilityBaseType(entropy_flux)
    , basis_functions_(basis_functions)
    , epsilon_(epsilon)
  {}

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& entity,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& eigenvectors,
                         const size_t dd) const override final
  {
    const VectorType slope = slope_limiter_.get(entity, stencil, eigenvectors, dd);
    const VectorType& u_bar = stencil[1];
    const FieldVector<VectorType, 2> reconstructed_values{u_bar - slope * 0.5, u_bar + slope * 0.5};

    VectorType thetas(0.);
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u = reconstructed_values[kk];
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        const auto& u0 = u[2 * ii];
        const auto& u1 = u[2 * ii + 1];
        const auto& ubar0 = u_bar[2 * ii];
        const auto& ubar1 = u_bar[2 * ii + 1];
        const auto& vj = basis_functions_.triangulation()[ii];
        const auto& vjplus1 = basis_functions_.triangulation()[ii + 1];
        FieldVector<RangeFieldType, 3> thetas_ii;
        if (!is_epsilon_realizable(ubar0, ubar1, vj, vjplus1, epsilon_)) {
          thetas[2 * ii] = 1.;
        } else {
          thetas_ii[0] = (epsilon_ - u0) / (ubar0 - u0);
          thetas_ii[1] =
              (u0 * vj - u1 + epsilon_ * std::sqrt(std::pow(vj, 2) + 1)) / ((ubar1 - u1) - (ubar0 - u0) * vj);
          thetas_ii[2] = (u0 * vjplus1 - u1 - epsilon_ * std::sqrt(std::pow(vjplus1, 2) + 1))
                         / ((ubar1 - u1) - (ubar0 - u0) * vjplus1);
          for (size_t ll = 0; ll < 3; ++ll)
            if (thetas_ii[ll] >= 0. && thetas_ii[ll] <= 1.)
              thetas[2 * ii] = std::max(thetas[2 * ii], thetas_ii[ll]);
        } // else (!realizable)
        thetas[2 * ii + 1] = thetas[2 * ii];
      } // ii
    } // local_reconstructed_values

    VectorType ret;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      assert(XT::Common::FloatCmp::le(thetas[ii], 1.) && XT::Common::FloatCmp::ge(thetas[ii], 0.));
      ret[ii] = slope[ii] * (1 - thetas[ii]);
    }
    this->ensure_solvability(entity, u_bar, ret);
    return ret;
  }

private:
  static bool is_epsilon_realizable(const RangeFieldType ubar0,
                                    const RangeFieldType ubar1,
                                    const RangeFieldType v0,
                                    const RangeFieldType v1,
                                    const RangeFieldType eps)
  {
    bool ret = (ubar0 >= eps) && (ubar1 <= v1 * ubar0 - eps * std::sqrt(std::pow(v1, 2) + 1))
               && (v0 * ubar0 + eps * std::sqrt(std::pow(v0, 2) + 1) <= ubar1);
    return ret;
  }

  const MomentBasis& basis_functions_;
  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
}; // class Dg1dRealizabilityLimitedSlope<...>


#if HAVE_QHULL
// Realizability limiter that ensures that the limited values are within the convex hull of the quadrature points. Uses
// single limiter variable for all components.
template <class GV,
          class MomentBasis,
          class EigenVectorWrapperType,
          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>>
class ConvexHullRealizabilityLimitedSlope
  : public SlopeBase<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType, 3>
  , public RealizabilityLimiterBase<GV, MomentBasis>
{
  using ThisType = ConvexHullRealizabilityLimitedSlope;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  static const size_t dimRange = MomentBasis::dimRange;
  using RealizabilityBaseType = RealizabilityLimiterBase<GV, MomentBasis>;
  using typename RealizabilityBaseType::E;
  using typename RealizabilityBaseType::EntropyFluxType;
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using typename BaseType::VectorType;
  using PlaneCoefficientsType = typename std::vector<std::pair<VectorType, RangeFieldType>>;

public:
  using typename BaseType::StencilType;

  ConvexHullRealizabilityLimitedSlope(const EntropyFluxType& entropy_flux,
                                      const MomentBasis& basis_functions,
                                      const RangeFieldType epsilon = 0.)
    : RealizabilityBaseType(entropy_flux)
    , basis_functions_(basis_functions)
    , epsilon_(epsilon)
  {
    calculate_plane_coefficients();
  }

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& entity,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& eigenvectors,
                         const size_t dd) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType& u_bar = stencil[1];
    if (!is_epsilon_realizable(u_bar, epsilon_))
      return zero_vector;
    const VectorType slope = slope_limiter_.get(entity, stencil, eigenvectors, dd);

    const FieldVector<VectorType, 2> reconstructed_values{u_bar - 0.5 * slope, u_bar + 0.5 * slope};

    RangeFieldType theta(-epsilon_);
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u = reconstructed_values[kk];
      // rescale u_l, u_bar
      auto u_bar_minus_u = u_bar - u;
      const auto factor = std::max(basis_functions_.density(u), basis_functions_.density(u_bar)) / 2.;
      u /= factor;
      u_bar_minus_u /= factor;

      for (const auto& coeffs : plane_coefficients_) {
        const auto& a = coeffs.first;
        const auto& b = coeffs.second;
        RangeFieldType theta_li = (b - a * u) / (a * u_bar_minus_u);
        if (XT::Common::FloatCmp::le(theta_li, 1.))
          theta = std::max(theta, theta_li);
      } // coeffs
    } // kk
    theta = std::min(epsilon_ + theta, 1.);

    assert(XT::Common::FloatCmp::le(theta, 1.) && XT::Common::FloatCmp::ge(theta, 0.));
    VectorType ret = slope * (1 - theta);
    this->ensure_solvability(entity, u_bar, ret);
    return ret;
  }

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
        plane_coefficients_[ii].first[jj] = *(facet.hyperplane().coordinates() + jj);
      plane_coefficients_[ii].second = -facet.hyperplane().offset();
    }
  } // void calculate_plane_coefficients()

  const MomentBasis& basis_functions_;
  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
  PlaneCoefficientsType plane_coefficients_;
}; // class ConvexHullRealizabilityLimitedSlope<...>


// Realizability limiter that ensures that the limited values are within the convex hull of the quadrature points. Uses
// single limiter variable for all components.
template <class GV,
          class MomentBasis,
          class EigenVectorWrapperType,

          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>>
class DgConvexHullRealizabilityLimitedSlope
  : public SlopeBase<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType, 3>
  , public RealizabilityLimiterBase<GV, MomentBasis>
{
  using ThisType = DgConvexHullRealizabilityLimitedSlope;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  static const size_t dimRange = MomentBasis::dimRange;
  static const size_t block_size = (MomentBasis::dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  using RealizabilityBaseType = RealizabilityLimiterBase<GV, MomentBasis>;
  using typename RealizabilityBaseType::E;
  using typename RealizabilityBaseType::EntropyFluxType;
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using typename BaseType::VectorType;
  using BlockRangeType = typename VectorType::BlockType;
  using BlockPlaneCoefficientsType = typename std::vector<std::pair<BlockRangeType, RangeFieldType>>;
  using PlaneCoefficientsType = FieldVector<BlockPlaneCoefficientsType, num_blocks>;

public:
  using typename BaseType::StencilType;

  DgConvexHullRealizabilityLimitedSlope(const EntropyFluxType& entropy_flux,
                                        const MomentBasis& basis_functions,
                                        const RangeFieldType epsilon = 0.)
    : RealizabilityBaseType(entropy_flux)
    , basis_functions_(basis_functions)
    , epsilon_(epsilon)
  {
    if (!basis_functions_.plane_coefficients()[0].size())
      basis_functions_.calculate_plane_coefficients();
    // replace the coefficients b by \tilde{b} = b - eps  \twonorm(a) to guarantee distance of eps to boundary
    plane_coefficients_ = basis_functions_.plane_coefficients();
    for (size_t jj = 0; jj < num_blocks; ++jj)
      for (auto& coeff : plane_coefficients_[jj])
        coeff.second -= epsilon_ * coeff.first.two_norm();
  }

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& entity,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& eigenvectors,
                         const size_t dd) const override final
  {
    const VectorType& u_bar = stencil[1];
    const VectorType slope = slope_limiter_.get(entity, stencil, eigenvectors, dd);
    const FieldVector<VectorType, 2> reconstructed_values{u_bar - slope * 0.5, u_bar + slope * 0.5};

    // vector to store thetas for each local reconstructed value
    FieldVector<RangeFieldType, num_blocks> thetas(0.);
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u = reconstructed_values[kk];
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        // Check realizability of u_bar in this block. The first condition avoids unnecessary repeated checking.
        if (thetas[jj] == 1. || !is_epsilon_realizable(u_bar.block(jj), jj)) {
          thetas[jj] = 1.;
          continue;
        }
        thetas[jj] = std::max(thetas[jj], get_block_theta(u.block(jj), u_bar.block(jj), jj));
      } // jj
    } // kk

    VectorType ret;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      assert(XT::Common::FloatCmp::le(thetas[jj], 1.) && XT::Common::FloatCmp::ge(thetas[jj], 0.));
      for (size_t ii = 0; ii < block_size; ++ii)
        ret.block(jj)[ii] = slope.block(jj)[ii] * (1 - thetas[jj]);
    }
    this->ensure_solvability(entity, u_bar, ret);
    return ret;
  } // ... get(...)

private:
  bool is_epsilon_realizable(const BlockRangeType& u, const size_t jj) const
  {
    for (const auto& coeff : plane_coefficients_[jj])
      if (u * coeff.first > coeff.second)
        return false;
    return true;
  }

  RangeFieldType get_block_theta(const BlockRangeType& u, const BlockRangeType& u_bar, const size_t jj) const
  {
    RangeFieldType theta(0.);
    auto u_bar_minus_u = u_bar - u;
    for (const auto& coeffs : plane_coefficients_[jj]) {
      const auto& a = coeffs.first;
      const auto& b = coeffs.second;
      RangeFieldType theta_li = (b - a * u) / (a * u_bar_minus_u);
      if (XT::Common::FloatCmp::le(theta_li, 1.))
        theta = std::max(theta, theta_li);
    } // coeffs
    return theta;
  } // ... get_block_theta(...)

  const MomentBasis& basis_functions_;
  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
  PlaneCoefficientsType plane_coefficients_;
}; // class DgConvexHullRealizabilityLimitedSlope<...>

#else // HAVE_QHULL

template <class GV,
          class MomentBasis,

          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>,
                                        FieldVector<typename MomentBasis::RangeFieldType, MomentBasis::dimRange>,
                                        EigenVectorWrapperType>>
class ConvexHullRealizabilityLimitedSlope
{
  static_assert(Dune::AlwaysFalse<MomentBasis>::value, "You are missing Qhull!");
};

template <class GV,
          class MomentBasis,

          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>,
                                        FieldVector<typename MomentBasis::RangeFieldType, MomentBasis::dimRange>,
                                        EigenVectorWrapperType>>
class DgConvexHullRealizabilityLimitedSlopeSlope
{
  static_assert(Dune::AlwaysFalse<MomentBasis>::value, "You are missing Qhull!");
};
#endif // HAVE_QHULL

#if HAVE_CLP
// Characteristic component-wise realizability limiter that ensures positivity of the components of u in
// noncharacteristic variables by solving a linear program.
template <class GV,
          class MomentBasis,
          class EigenVectorWrapperType,
          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>>
class LpPositivityLimitedSlope
  : public SlopeBase<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType, 3>
  , public RealizabilityLimiterBase<GV, MomentBasis>
{
  using ThisType = LpPositivityLimitedSlope;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  static const size_t dimRange = MomentBasis::dimRange;
  using RealizabilityBaseType = RealizabilityLimiterBase<GV, MomentBasis>;
  using typename RealizabilityBaseType::E;
  using typename RealizabilityBaseType::EntropyFluxType;
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;

public:
  using typename BaseType::StencilType;

  LpPositivityLimitedSlope(const EntropyFluxType& entropy_flux, const RangeFieldType epsilon)
    : RealizabilityBaseType(entropy_flux)
    , epsilon_(epsilon)
  {}

  LpPositivityLimitedSlope(const ThisType& other)
    : BaseType(other)
    , RealizabilityBaseType(other)
    , epsilon_(other.epsilon_)
    , A_tilde_transposed_(std::make_unique<MatrixType>(dimRange, dimRange))
    , slope_limiter_(other.slope_limiter_)
  {}

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& entity,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& eigenvectors,
                         const size_t dd) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType& u_bar = stencil[1];
    if (!is_epsilon_realizable(u_bar, epsilon_))
      return zero_vector;
    const VectorType slope_char = slope_limiter_.get_char(entity, stencil, eigenvectors, dd);
    if (XT::Common::FloatCmp::eq(slope_char, zero_vector))
      return zero_vector;
    FieldVector<VectorType, 2> thetas;
    VectorType u_bar_char;
    eigenvectors.apply_inverse_eigenvectors(dd, u_bar, u_bar_char);
    const FieldVector<VectorType, 2> reconstructed_values{u_bar_char - 0.5 * slope_char, u_bar_char + 0.5 * slope_char};
    auto& A_tilde_transposed = *A_tilde_transposed_;
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u_char = reconstructed_values[kk];
      thetas[kk] = solve_linear_program(u_char, u_bar_char, eigenvectors.eigenvectors(dd), A_tilde_transposed);
    } // kk

    VectorType slope;
    eigenvectors.apply_eigenvectors(dd, slope_char, slope);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      const auto theta_max_ii = std::max(thetas[0][ii], thetas[1][ii]);
      slope[ii] = slope[ii] * (1 - theta_max_ii);
    }
    this->ensure_solvability(entity, u_bar, slope);
    return slope;
  }

private:
  static bool is_epsilon_realizable(const VectorType& u, const RangeFieldType epsilon)
  {
    for (const auto& val : u)
      if (val < epsilon)
        return false;
    return true;
  }

  VectorType solve_linear_program(const VectorType& u_char,
                                  const VectorType u_bar_char,
                                  const MatrixType& A,
                                  MatrixType& A_tilde_transposed) const
  {
    // calculate data
    VectorType u;
    A.mv(u_char, u);
    VectorType u_minus_u_bar = u_char - u_bar_char;
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        A_tilde_transposed.set_entry(jj, ii, A.get_entry(ii, jj) * u_minus_u_bar[jj]);

    // Create LP with dimRange rows and columns */
    constexpr int num_rows = static_cast<int>(dimRange);
    constexpr int num_cols = static_cast<int>(dimRange);
    auto lp_ptr = std::make_unique<ClpSimplex>(false);
    auto& lp = *lp_ptr;

    // set number of rows, columns will be added below
    lp.resize(num_rows, 0);

    // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
    // 0..num_rows
    std::array<int, num_rows> row_indices;
    for (int ii = 0; ii < num_rows; ++ii)
      row_indices[ii] = ii;

    // set columns
    for (int jj = 0; jj < num_cols; ++jj)
      lp.addColumn(num_rows, row_indices.data(), &(A_tilde_transposed[jj][0]), 0., 1., 1.);

    lp.setLogLevel(0);

    lp.setMaximumWallSeconds(60);

    // set rhs
    for (int ii = 0; ii < num_rows; ++ii)
      lp.setRowUpper(ii, u[ii] - epsilon_);

    // Now solve
    VectorType ret;
    lp.primal();
    const auto* thetas = lp.primalColumnSolution();
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = thetas[ii];
    if (!lp.isProvenOptimal())
      std::fill(ret.begin(), ret.end(), 1.);
    return ret;
  } // void solve_linear_program(...)

  const RangeFieldType epsilon_;
  mutable std::unique_ptr<MatrixType> A_tilde_transposed_;
  const SlopeType slope_limiter_;
}; // class LpPositivityLimitedSlope<...>

template <class GV, class MomentBasis, class EigenVectorWrapperType, class SlopeType>
constexpr size_t LpPositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType, SlopeType>::dimRange;


// Realizability limiter that solves a linear program to ensure the reconstructed values are still in the numerically
// realizable set, i.e. in the convex hull of basis evaluations.
template <class GV,
          class MomentBasis,
          class EigenVectorWrapperType,
          class SlopeType = MinmodSlope<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>>
class LpConvexhullRealizabilityLimitedSlope
  : public SlopeBase<XT::Grid::extract_entity_t<GV>, EigenVectorWrapperType>
  , public RealizabilityLimiterBase<GV, MomentBasis>
{
  using ThisType = LpConvexhullRealizabilityLimitedSlope;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  static_assert(dimRange < std::numeric_limits<int>::max(), "");
  static constexpr size_t num_rows = dimRange;
  using RealizabilityBaseType = RealizabilityLimiterBase<GV, MomentBasis>;
  using typename RealizabilityBaseType::E;
  using typename RealizabilityBaseType::EntropyFluxType;
  using BaseType = SlopeBase<E, EigenVectorWrapperType, 3>;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;

public:
  using typename BaseType::StencilType;

  LpConvexhullRealizabilityLimitedSlope(const EntropyFluxType& entropy_flux,
                                        const MomentBasis& basis_functions,
                                        const RangeFieldType epsilon)
    : RealizabilityBaseType(entropy_flux)
    , epsilon_(epsilon)
    , basis_functions_(basis_functions)
    , basis_values_(
          std::make_shared<std::vector<VectorType>>(XT::Data::merged_quadrature(basis_functions_.quadratures()).size()))
  {
    const auto& quadrature = XT::Data::merged_quadrature(basis_functions_.quadratures());
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      (*basis_values_)[ii] = basis_functions.evaluate(quadrature[ii].position());
  }

  LpConvexhullRealizabilityLimitedSlope(const LpConvexhullRealizabilityLimitedSlope& other)
    : BaseType(other)
    , RealizabilityBaseType(other)
    , epsilon_(other.epsilon_)
    , basis_functions_(other.basis_functions_)
    , basis_values_(other.basis_values_)
    , lp_(nullptr)
    , A_tilde_transposed_(std::make_unique<MatrixType>(dimRange, dimRange))
    , slope_limiter_(other.slope_limiter_)
  {}

  virtual BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  virtual VectorType get(const E& entity,
                         const StencilType& stencil,
                         const EigenVectorWrapperType& eigenvectors,
                         const size_t dd) const override final
  {
    VectorType slope_char = slope_limiter_.get_char(entity, stencil, eigenvectors, dd);
    static const VectorType zero_vector(0.);
    if (XT::Common::FloatCmp::eq(slope_char, zero_vector))
      return zero_vector;
    FieldVector<VectorType, 2> thetas;
    const VectorType& u_bar = stencil[1];
    VectorType u_bar_char;
    eigenvectors.apply_inverse_eigenvectors(dd, u_bar, u_bar_char);
    const FieldVector<VectorType, 2> reconstructed_values_char{u_bar_char - 0.5 * slope_char,
                                                               u_bar_char + 0.5 * slope_char};
    auto& A_tilde_transposed = *A_tilde_transposed_;
    for (size_t kk = 0; kk < reconstructed_values_char.size(); ++kk)
      thetas[kk] = solve_linear_program(
          reconstructed_values_char[kk], u_bar_char, eigenvectors.eigenvectors(dd), A_tilde_transposed);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      auto theta_max_ii = std::max(thetas[0][ii], thetas[1][ii]);
      if (theta_max_ii > 0.)
        theta_max_ii = std::min(1., theta_max_ii + epsilon_);
      slope_char[ii] = slope_char[ii] * (1 - theta_max_ii);
    }
    // Ensure positive density
    // For that purpose, get slope in ordinary coordinates
    VectorType slope;
    eigenvectors.apply_eigenvectors(dd, slope_char, slope);
    const FieldVector<VectorType, 2> reconstructed_values{u_bar - 0.5 * slope, u_bar + 0.5 * slope};
    RangeFieldType theta_pos(0);
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const auto& u = reconstructed_values[kk];
      const auto density_u = basis_functions_.density(u);
      const auto density_u_bar = basis_functions_.density(u_bar);
      if (density_u >= density_u_bar)
        continue;
      else if (density_u < epsilon_)
        theta_pos = std::max(theta_pos, (epsilon_ - density_u) / (density_u_bar - density_u));
    } // kk
    slope *= 1 - theta_pos;
    this->ensure_solvability(entity, u_bar, slope);
    return slope;
  }

private:
  void setup_linear_program() const
  {
    // We start with creating a model with dimRange rows and num_quad_points+1 columns */
    assert(basis_values_->size() + dimRange < std::numeric_limits<int>::max());
    size_t num_cols = basis_values_->size() + dimRange; /* variables are x_1, ..., x_{num_quad_points}, theta_1,
                                                                            ..., theta_{dimRange} */
    lp_ = std::make_unique<ClpSimplex>(false);
    auto& lp = *lp_;
    // set number of rows
    lp.resize(num_rows, 0);

    // set columns for quadrature points
    assert(basis_values_->size() == num_cols - dimRange);
    static auto row_indices = create_row_indices();
    for (size_t ii = 0; ii < num_cols - dimRange; ++ii) {
      // First argument: number of elements in column
      // Second/Third argument: indices/values of column entries
      // Fourth/Fifth argument: lower/upper column bound, i.e. lower/upper bound for x_i. As all x_i should be
      // positive, set to 0/inf, which is the default.
      // Sixth argument: Prefactor in objective for x_i, this is 0 for all x_i, which is also the default;
      lp.addColumn(static_cast<int>(num_rows), row_indices.data(), &((*basis_values_)[ii][0]));
    }

    // add theta columns (set to random values, will be set correctly in solve_linear_program)
    // The bounds for theta should be [0,1]. Also sets the prefactor in the objective to 1 for the thetas.
    for (size_t ii = 0; ii < dimRange; ++ii)
      lp.addColumn(static_cast<int>(num_rows), row_indices.data(), &((*basis_values_)[0][0]), 0., 1., 1.);
    lp.setLogLevel(0);
  } // void setup_linear_program()

  VectorType solve_linear_program(const VectorType& u_char,
                                  const VectorType& u_bar_char,
                                  const MatrixType& A,
                                  MatrixType& A_tilde_transposed) const
  {
    // calculate data
    VectorType u_minus_u_bar_char = u_char - u_bar_char;
    VectorType u;
    A.mv(u_char, u);
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        A_tilde_transposed.set_entry(jj, ii, A.get_entry(ii, jj) * u_minus_u_bar_char[jj]);

    // setup linear program
    setup_linear_program();
    auto& lp = *lp_;
    size_t num_cols = basis_values_->size() + dimRange; // see above

    // set rhs (equality constraints, so set both bounds equal)
    for (size_t ii = 0; ii < num_rows; ++ii) {
      lp.setRowLower(static_cast<int>(ii), u[ii]);
      lp.setRowUpper(static_cast<int>(ii), u[ii]);
    }

    // delete old theta columns.
    FieldVector<int, dimRange> theta_columns;
    for (size_t ii = 0; ii < dimRange; ++ii)
      theta_columns[ii] = static_cast<int>(num_cols - dimRange + ii);
    lp.deleteColumns(dimRange, &(theta_columns[0]));

    // set new theta columns
    static auto row_indices = create_row_indices();
    for (size_t jj = 0; jj < dimRange; ++jj)
      lp.addColumn(num_rows, row_indices.data(), &(A_tilde_transposed[jj][0]), 0., 1., 1.);

    // Now solve
    lp.setMaximumWallSeconds(60);
    lp.primal();
    const auto* thetas_ptr = lp.primalColumnSolution();
    VectorType thetas;
    for (size_t ii = 0; ii < dimRange; ++ii)
      thetas[ii] = thetas_ptr[ii];
    if (!lp.isProvenOptimal())
      std::fill(thetas.begin(), thetas.end(), 1.);
    return thetas;
  }

  // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
  // 0..num_rows
  static std::array<int, num_rows> create_row_indices()
  {
    std::array<int, dimRange> ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = static_cast<int>(ii);
    return ret;
  }

  const RangeFieldType epsilon_;
  const MomentBasis& basis_functions_;
  std::shared_ptr<std::vector<VectorType>> basis_values_;
  mutable std::unique_ptr<ClpSimplex> lp_;
  mutable std::unique_ptr<MatrixType> A_tilde_transposed_;
  const SlopeType slope_limiter_;
};

template <class GV, class MomentBasis, class EigenVectorWrapperType, class SlopeType>
constexpr size_t LpConvexhullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType, SlopeType>::dimRange;

#endif // HAVE_CLP


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_SLOPES_HH
