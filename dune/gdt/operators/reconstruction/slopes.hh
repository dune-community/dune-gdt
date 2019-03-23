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

#include "internal.hh"

namespace Dune {
namespace GDT {


template <class VectorType, class MatrixType, size_t stencil_size = 3>
class SlopeBase
{
public:
  virtual ~SlopeBase() = default;

  using StencilType = FieldVector<VectorType, stencil_size>;

  // returns (limited) calculated slope in characteristic coordinates
  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& eigenvectors) const = 0;
};


// Central slope without any limiting
template <class VectorType, class MatrixType, size_t stencil_size>
class CentralSlope : public SlopeBase<VectorType, MatrixType, stencil_size>
{
  using BaseType = SlopeBase<VectorType, MatrixType, stencil_size>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& stencil_char,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    const VectorType& u_left = stencil_char[0];
    const VectorType& u_right = stencil_char[2];
    return (u_right - u_left) / 2.;
  }
};


// Left slope without any limiting
template <class VectorType, class MatrixType, size_t stencil_size>
class LeftSlope : public SlopeBase<VectorType, MatrixType, stencil_size>
{
  using BaseType = SlopeBase<VectorType, MatrixType, stencil_size>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& stencil_char,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    const VectorType& u_left = stencil_char[0];
    const VectorType& u = stencil_char[1];
    return u - u_left;
  }
};


// Right slope without any limiting
template <class VectorType, class MatrixType, size_t stencil_size>
class RightSlope : public SlopeBase<VectorType, MatrixType, stencil_size>
{
  using BaseType = SlopeBase<VectorType, MatrixType, stencil_size>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& stencil_char,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    const VectorType& u = stencil_char[1];
    const VectorType& u_right = stencil_char[2];
    return u_right - u;
  }
};


// Zero slope
template <class VectorType, class MatrixType, size_t stencil_size>
class NoSlope : public SlopeBase<VectorType, MatrixType, stencil_size>
{
  using BaseType = SlopeBase<VectorType, MatrixType, stencil_size>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& /*stencil_char*/,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    return VectorType(0.);
  }
};


template <class VectorType, class MatrixType>
class MinmodSlope : public SlopeBase<VectorType, MatrixType, 3>
{
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& stencil_char,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    const VectorType slope_left = stencil_char[1] - stencil_char[0];
    const VectorType slope_right = stencil_char[2] - stencil_char[1];
    return minmod(slope_left, slope_right);
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


template <class VectorType, class MatrixType>
class McSlope : public SlopeBase<VectorType, MatrixType, 3>
{
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;
  using MinmodType = MinmodSlope<VectorType, MatrixType>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& stencil_char,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    const VectorType& u_left = stencil_char[0];
    const VectorType& u = stencil_char[1];
    const VectorType& u_right = stencil_char[2];
    const VectorType slope_left_twice = (u - u_left) * 2.;
    const VectorType slope_right_twice = (u_right - u) * 2.;
    const VectorType slope_center = (u_right - u_left) / 2.;
    const VectorType first_slope = MinmodType::minmod(slope_left_twice, slope_right_twice);
    return MinmodType::minmod(first_slope, slope_center);
  }
};


template <class VectorType, class MatrixType>
class SuperbeeSlope : public SlopeBase<VectorType, MatrixType, 3>
{
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;
  using MinmodType = MinmodSlope<VectorType, MatrixType>;

public:
  using typename BaseType::StencilType;

  virtual VectorType get(const StencilType& /*stencil*/,
                         const StencilType& stencil_char,
                         const MatrixType& /*eigenvectors*/) const override final
  {
    const VectorType slope_left = stencil_char[1] - stencil_char[0];
    const VectorType slope_right = stencil_char[2] - stencil_char[1];
    const VectorType first_slope = MinmodType::minmod(slope_left, slope_right * 2.);
    const VectorType second_slope = MinmodType::minmod(slope_left * 2., slope_right);
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


// Realizability limiter that ensures positivity of the components of u in noncharacteristic variables. Uses single
// limiter variable for all components.
// TODO: Make usable with interface quadratures different from the midpoint quadrature
// See dune/gdt/operators/fv/entropybased/realizability.hh
template <class RangeFieldType,
          size_t dimRange,
          class MatrixType,
          class SlopeType = MinmodSlope<FieldVector<RangeFieldType, dimRange>, MatrixType>>
class PositivityLimitedSlope : public SlopeBase<FieldVector<RangeFieldType, dimRange>, MatrixType, 3>
{
  using VectorType = FieldVector<RangeFieldType, dimRange>;
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;

public:
  using typename BaseType::StencilType;

  // This limiter ensures u_i >= epsilon for all components u_i of u.
  PositivityLimitedSlope(const RangeFieldType epsilon = 0.)
    : epsilon_(epsilon)
  {}

  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& A) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType& u_bar = stencil[1];
    if (!is_epsilon_realizable(u_bar, epsilon_))
      return zero_vector;
    const VectorType slope = slope_limiter_.get(stencil, stencil_char, A);
    // this needs to be changed for other interface quadratures (see above)
    const VectorType& u_bar_char = stencil_char[1];
    const FieldVector<VectorType, 2> reconstructed_values{u_bar_char - 0.5 * slope, u_bar_char + 0.5 * slope};
    VectorType thetas(0.);
    VectorType u;
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u_char = reconstructed_values[kk];
      // convert back to ordinary coordinates
      A.mv(u_char, u);
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
// TODO: Make usable with interface quadratures different from the midpoint quadrature
// See dune/gdt/operators/fv/entropybased/realizability.hh
template <class RangeFieldType,
          size_t dimRange,
          class MatrixType,
          class SlopeType = MinmodSlope<XT::Common::BlockedFieldVector<RangeFieldType, dimRange / 2, 2>, MatrixType>>
class Dg1dRealizabilityLimitedSlope
  : public SlopeBase<XT::Common::BlockedFieldVector<RangeFieldType, dimRange / 2, 2>, MatrixType, 3>
{
  using VectorType = XT::Common::BlockedFieldVector<RangeFieldType, dimRange / 2, 2>;
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;
  using BasisfunctionType = Dune::GDT::PartialMomentBasis<RangeFieldType, 1, RangeFieldType, dimRange, 1, 1, 1>;

public:
  using typename BaseType::StencilType;

  // This limiter ensures u_i >= epsilon for all components u_i of u.
  Dg1dRealizabilityLimitedSlope(const BasisfunctionType& basis_functions, const RangeFieldType epsilon = 0.)
    : basis_functions_(basis_functions)
    , epsilon_(epsilon)
  {}

  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& A) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType slope = slope_limiter_.get(stencil, stencil_char, A);
    // this needs to be changed for other interface quadratures (see above)
    const VectorType& u_bar_char = stencil_char[1];
    const FieldVector<VectorType, 2> reconstructed_values{u_bar_char - slope * 0.5, u_bar_char + slope * 0.5};
    const VectorType& u_bar = stencil[1];

    VectorType thetas(0.);
    VectorType u;
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u_char = reconstructed_values[kk];
      // convert back to ordinary coordinates
      A.mv(u_char, u);
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

  const BasisfunctionType& basis_functions_;
  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
}; // class Dg1dRealizabilityLimitedSlope<...>


#if HAVE_QHULL
// Realizability limiter that ensures that the limited values are within the convex hull of the quadrature points. Uses
// single limiter variable for all components.
// TODO: Make usable with interface quadratures different from the midpoint quadrature
// See dune/gdt/operators/fv/entropybased/realizability.hh
template <
    class BasisfunctionType,
    class MatrixType,
    class SlopeType =
        MinmodSlope<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>, MatrixType>>
class ConvexHullRealizabilityLimitedSlope
  : public SlopeBase<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>,
                     MatrixType,
                     3>
{
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  static const size_t dimRange = BasisfunctionType::dimRange;
  using VectorType = FieldVector<RangeFieldType, dimRange>;
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;
  using PlaneCoefficientsType = typename std::vector<std::pair<VectorType, RangeFieldType>>;

public:
  using typename BaseType::StencilType;

  ConvexHullRealizabilityLimitedSlope(const BasisfunctionType& basis_functions, const RangeFieldType epsilon = 0.)
    : basis_functions_(basis_functions)
    , epsilon_(epsilon)
  {
    calculate_plane_coefficients();
  }

  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& A) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType& u_bar = stencil[1];
    if (!is_epsilon_realizable(u_bar, epsilon_))
      return zero_vector;
    const VectorType slope = slope_limiter_.get(stencil, stencil_char, A);

    // this needs to be changed for other interface quadratures (see above)
    const VectorType& u_bar_char = stencil_char[1];
    const FieldVector<VectorType, 2> reconstructed_values{u_bar_char - 0.5 * slope, u_bar_char + 0.5 * slope};

    RangeFieldType theta(-epsilon_);
    VectorType u;
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u_char = reconstructed_values[kk];
      // convert back to ordinary coordinates
      A.mv(u_char, u);
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

  const BasisfunctionType& basis_functions_;
  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
  PlaneCoefficientsType plane_coefficients_;
}; // class ConvexHullRealizabilityLimitedSlope<...>


// Realizability limiter that ensures that the limited values are within the convex hull of the quadrature points. Uses
// single limiter variable for all components.
// TODO: Make usable with interface quadratures different from the midpoint quadrature
// See dune/gdt/operators/fv/entropybased/realizability.hh
template <
    class BasisfunctionType,
    class MatrixType,
    class SlopeType =
        MinmodSlope<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>, MatrixType>>
class DgConvexHullRealizabilityLimitedSlope
  : public SlopeBase<
        XT::Common::BlockedFieldVector<typename BasisfunctionType::RangeFieldType,
                                       BasisfunctionType::dimRange / ((BasisfunctionType::dimDomain == 1) ? 2 : 4),
                                       (BasisfunctionType::dimDomain == 1) ? 2 : 4>,
        MatrixType,
        3>
{
  using ThisType = DgConvexHullRealizabilityLimitedSlope;
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  static const size_t dimRange = BasisfunctionType::dimRange;
  static const size_t block_size = (BasisfunctionType::dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  using VectorType = XT::Common::BlockedFieldVector<RangeFieldType, num_blocks, block_size>;
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;
  using BlockRangeType = typename VectorType::BlockType;
  using BlockPlaneCoefficientsType = typename std::vector<std::pair<BlockRangeType, RangeFieldType>>;
  using PlaneCoefficientsType = FieldVector<BlockPlaneCoefficientsType, num_blocks>;

public:
  using typename BaseType::StencilType;

  DgConvexHullRealizabilityLimitedSlope(const BasisfunctionType& basis_functions, const RangeFieldType epsilon = 0.)
    : basis_functions_(basis_functions)
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

  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& A) const override final
  {
    const VectorType& u_bar = stencil[1];
    const VectorType slope = slope_limiter_.get(stencil, stencil_char, A);

    // this needs to be changed for other interface quadratures (see above)
    const VectorType& u_bar_char = stencil_char[1];
    const FieldVector<VectorType, 2> reconstructed_values{u_bar_char - slope * 0.5, u_bar_char + slope * 0.5};

    // vector to store thetas for each local reconstructed value
    FieldVector<RangeFieldType, num_blocks> thetas(0.);
    VectorType u;
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u_char = reconstructed_values[kk];
      // convert back to ordinary coordinates
      A.mv(u_char, u);
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

  const BasisfunctionType& basis_functions_;
  const RangeFieldType epsilon_;
  const SlopeType slope_limiter_;
  PlaneCoefficientsType plane_coefficients_;
}; // class DgConvexHullRealizabilityLimitedSlope<...>

#else // HAVE_QHULL

template <
    class BasisfunctionType,
    class MatrixType,
    class SlopeType =
        MinmodSlope<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>, MatrixType>>
class ConvexHullRealizabilityLimitedSlope
{
  static_assert(Dune::AlwaysFalse<BasisfunctionType>::value, "You are missing Qhull!");
};

template <
    class BasisfunctionType,
    class MatrixType,
    class SlopeType =
        MinmodSlope<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>, MatrixType>>
class DgConvexHullRealizabilityLimitedSlopeSlope
{
  static_assert(Dune::AlwaysFalse<BasisfunctionType>::value, "You are missing Qhull!");
};
#endif // HAVE_QHULL

#if HAVE_CLP
// Characteristic component-wise realizability limiter that ensures positivity of the components of u in
// noncharacteristic variables by solving a linear program.
// TODO: Make usable with interface quadratures different from the midpoint quadrature
// See dune/gdt/operators/fv/entropybased/realizability.hh
template <class RangeFieldType,
          size_t dimRange,
          class MatrixType,
          class SlopeType = MinmodSlope<FieldVector<RangeFieldType, dimRange>, MatrixType>>
class LpPositivityLimitedSlope : public SlopeBase<FieldVector<RangeFieldType, dimRange>, MatrixType, 3>
{
  using VectorType = FieldVector<RangeFieldType, dimRange>;
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;

public:
  using typename BaseType::StencilType;

  LpPositivityLimitedSlope(const RangeFieldType epsilon)
    : epsilon_(epsilon)
  {}

  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& A) const override final
  {
    static const VectorType zero_vector(0.);
    const VectorType& u_bar = stencil[1];
    if (!is_epsilon_realizable(u_bar, epsilon_))
      return zero_vector;
    const VectorType slope = slope_limiter_.get(stencil, stencil_char, A);
    if (XT::Common::FloatCmp::eq(slope, zero_vector))
      return zero_vector;
    FieldVector<VectorType, 2> thetas;
    // this needs to be changed for other interface quadratures (see above)
    const VectorType& u_bar_char = stencil_char[1];
    const FieldVector<VectorType, 2> reconstructed_values{u_bar_char - 0.5 * slope, u_bar_char + 0.5 * slope};
    auto& A_tilde_transposed = *A_tilde_transposed_;
    for (size_t kk = 0; kk < reconstructed_values.size(); ++kk) {
      const VectorType& u_char = reconstructed_values[kk];
      thetas[kk] = solve_linear_program(u_char, u_bar_char, A, A_tilde_transposed);
    } // kk

    VectorType ret;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      const auto theta_max_ii = std::max(thetas[0][ii], thetas[1][ii]);
      ret[ii] = slope[ii] * (1 - theta_max_ii);
    }
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
        A_tilde_transposed[jj][ii] = A[ii][jj] * u_minus_u_bar[jj];

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
  mutable XT::Common::PerThreadValue<MatrixType> A_tilde_transposed_;
  const SlopeType slope_limiter_;
}; // class LpPositivityLimitedSlope<...>


// Realizability limiter that solves a linear program to ensure the reconstructed values are still in the numerically
// realizable set, i.e. in the convex hull of basis evaluations.
// TODO: Make usable with interface quadratures different from the midpoint quadrature
// See dune/gdt/operators/fv/entropybased/realizability.hh
template <
    class BasisfunctionType,
    class MatrixType,
    class SlopeType =
        MinmodSlope<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>, MatrixType>>
class LpConvexhullRealizabilityLimitedSlope
  : public SlopeBase<FieldVector<typename BasisfunctionType::RangeFieldType, BasisfunctionType::dimRange>, MatrixType>
{
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  static constexpr size_t dimRange = BasisfunctionType::dimRange;
  static_assert(dimRange < std::numeric_limits<int>::max(), "");
  static constexpr size_t num_rows = dimRange;
  using VectorType = FieldVector<RangeFieldType, dimRange>;
  using BaseType = SlopeBase<VectorType, MatrixType, 3>;

public:
  using typename BaseType::StencilType;

  LpConvexhullRealizabilityLimitedSlope(const BasisfunctionType& basis_functions, const RangeFieldType epsilon)
    : epsilon_(epsilon)
    , basis_functions_(basis_functions)
    , basis_values_(basis_functions_.quadratures().merged().size())
  {
    const auto& quadrature = basis_functions.quadratures().merged();
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      basis_values_[ii] = basis_functions.evaluate(quadrature[ii].position());
  }

  virtual VectorType
  get(const StencilType& stencil, const StencilType& stencil_char, const MatrixType& A) const override final
  {
    const VectorType slope_char = slope_limiter_.get(stencil, stencil_char, A);
    static const VectorType zero_vector(0.);
    if (XT::Common::FloatCmp::eq(slope_char, zero_vector))
      return zero_vector;
    FieldVector<VectorType, 2> thetas;
    // this needs to be changed for other interface quadratures (see above)
    const auto& u_bar_char = stencil_char[1];
    const FieldVector<VectorType, 2> reconstructed_values_char{u_bar_char - 0.5 * slope_char,
                                                               u_bar_char + 0.5 * slope_char};
    auto& A_tilde_transposed = *A_tilde_transposed_;
    for (size_t kk = 0; kk < reconstructed_values_char.size(); ++kk)
      thetas[kk] = solve_linear_program(reconstructed_values_char[kk], u_bar_char, A, A_tilde_transposed);
    VectorType ret;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      auto theta_max_ii = std::max(thetas[0][ii], thetas[1][ii]);
      if (theta_max_ii > 0.)
        theta_max_ii = std::min(1., theta_max_ii + epsilon_);
      ret[ii] = slope_char[ii] * (1 - theta_max_ii);
    }
    // Ensure positive density
    // For that purpose, get slope in ordinary coordinates
    VectorType slope;
    A.mv(ret, slope);
    const auto& u_bar = stencil[1];
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
    ret *= 1 - theta_pos;
    return ret;
  }

private:
  void setup_linear_program() const
  {
    // It should be possible to reuse the same ClpSimplex class over and over again, only changing the relevant columns.
    // However, in rare cases the class seems to get corrupted, always claiming the LP is infeasible, even if it is not.
    // Creating a new LP from time to time seems to fix this problem.
    thread_local int counter;
    ++counter;
    if (!*lp_ || !(counter % 100)) {
      counter = 0;
      // We start with creating a model with dimRange rows and num_quad_points+1 columns */
      assert(basis_values_.size() + dimRange < std::numeric_limits<int>::max());
      size_t num_cols = basis_values_.size() + dimRange; /* variables are x_1, ..., x_{num_quad_points}, theta_1,
                                                                              ..., theta_{dimRange} */
      *lp_ = std::make_unique<ClpSimplex>(false);
      auto& lp = **lp_;
      // set number of rows
      lp.resize(num_rows, 0);

      // set columns for quadrature points
      assert(basis_values_.size() == num_cols - dimRange);
      static auto row_indices = create_row_indices();
      for (size_t ii = 0; ii < num_cols - dimRange; ++ii) {
        // First argument: number of elements in column
        // Second/Third argument: indices/values of column entries
        // Fourth/Fifth argument: lower/upper column bound, i.e. lower/upper bound for x_i. As all x_i should be
        // positive, set to 0/inf, which is the default.
        // Sixth argument: Prefactor in objective for x_i, this is 0 for all x_i, which is also the default;
        lp.addColumn(static_cast<int>(num_rows), row_indices.data(), &(basis_values_[ii][0]));
      }

      // add theta columns (set to random values, will be set correctly in solve_linear_program)
      // The bounds for theta should be [0,1]. Also sets the prefactor in the objective to 1 for the thetas.
      for (size_t ii = 0; ii < dimRange; ++ii)
        lp.addColumn(static_cast<int>(num_rows), row_indices.data(), &(basis_values_[0][0]), 0., 1., 1.);
      lp.setLogLevel(0);
    } // if (!lp_)
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
        A_tilde_transposed[jj][ii] = A[ii][jj] * u_minus_u_bar_char[jj];

    // setup linear program
    setup_linear_program();
    auto& lp = **lp_;
    size_t num_cols = basis_values_.size() + dimRange; // see above

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
  const BasisfunctionType& basis_functions_;
  std::vector<VectorType> basis_values_;
  mutable XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>> lp_;
  mutable XT::Common::PerThreadValue<MatrixType> A_tilde_transposed_;
  const SlopeType slope_limiter_;
};
#endif // HAVE_CLP


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_SLOPES_HH
