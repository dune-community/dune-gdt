// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_SPHERICALHARMONICS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_SPHERICALHARMONICS_HH

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <dune/xt/common/coordinates.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


// TODO: use complex arithmetic, currently only usable for Pn Models in 2D, test for only_positive = false
template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_positive = true>
class SphericalHarmonics
    : public BasisfunctionsInterface<DomainFieldType,
                                     3,
                                     RangeFieldType,
                                     only_positive ? ((order + 1) * (order + 2)) / 2 : (order + 1) * (order + 1),
                                     1,
                                     fluxDim>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimRange = only_positive ? ((order + 1) * (order + 2)) / 2 : (order + 1) * (order + 1);
  static const size_t dimFlux = fluxDim;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, dimFlux> BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  virtual RangeType evaluate(const DomainType& v) const override
  {
    const auto v_spherical = XT::Common::CoordinateConverter<DomainFieldType>::to_spherical(v);
    return evaluate_in_spherical_coords(v_spherical);
  } // ... evaluate(...)

  RangeType evaluate_in_spherical_coords(const FieldVector<DomainFieldType, 2>& coords) const
  {
    const DomainFieldType theta = coords[0];
    const DomainFieldType phi = coords[1];
    RangeType ret(0);
    // TODO: use complex arithmetic, remove real() call
    assert(order <= std::numeric_limits<int>::max());
    for (unsigned int ll = 0; ll <= static_cast<unsigned int>(order); ++ll)
      for (int mm = only_positive ? 0 : -static_cast<int>(ll); mm <= static_cast<int>(ll); ++mm)
        ret[helper<only_positive>::pos(ll, mm)] = boost::math::spherical_harmonic(ll, mm, theta, phi).real();
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override
  {
    RangeType ret(0);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 1;
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    return mass_matrix();
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimFlux> ret(MatrixType(dimRange, dimRange, 0));
    ret[0] = create_Bx();
    ret[1] = create_Bz();
    //    if (dimFlux == 3)
    //      ret[2] = create_By();
    return ret;
  } // ... mass_matrix_with_v()

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      component_visualizer<DiscreteFunctionType, dimRange, 0>(u_n, filename_prefix, ii, std::sqrt(4 * M_PI));
    };
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) { return XT::Common::to_string(val[0] * std::sqrt(4 * M_PI), 15); };
  } // ... stringifier()

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeType u_iso(0), alpha_iso(0);
    u_iso[0] = u[0];
    alpha_iso[0] = std::log(u[0] / (4. * M_PI));
    return std::make_pair(u_iso, alpha_iso);
  }

private:
  static RangeFieldType A_lm(const size_t l, const int m)
  {
    return std::sqrt((l + m) * (l - m) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static RangeFieldType B_lm(const size_t l, const int m)
  {
    return std::sqrt((l + m) * (l + m - 1.) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static MatrixType create_Bx()
  {
    MatrixType Bx(dimRange, dimRange, 0);
    const auto& pos = helper<only_positive>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = only_positive ? 0 : -int(l1); size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -int(l2); size_t(std::abs(m2)) <= l2; ++m2) {
            size_t row = pos(l1, m1);
            size_t col = pos(l2, only_positive ? std::abs(m2) : m2);
            RangeFieldType factor = !only_positive ? 1. : (m2 < 0 ? std::pow(-1., m2) : 1.);
            if (l1 == l2 + 1 && m1 == m2 + 1)
              Bx[row][col] += -0.5 * factor * B_lm(l2 + 1, m2 + 1);
            if (l1 == l2 - 1 && m1 == m2 + 1)
              Bx[row][col] += 0.5 * factor * B_lm(l2, -m2);
            if (l1 == l2 + 1 && m1 == m2 - 1)
              Bx[row][col] += 0.5 * factor * B_lm(l2 + 1, -m2 - 1);
            if (l1 == l2 - 1 && m1 == m2 - 1)
              Bx[row][col] += -0.5 * factor * B_lm(l2, m2);
          } // m2
        } // l2
      } // m1
    } // l1
    return Bx;
  } // ... create_Bx()

  //    static MatrixType create_By()
  //    {
  //      MatrixType By(dimRange, dimRange, 0);
  //      const auto& pos = helper<only_positive>::pos;
  //      for (size_t l1 = 0; l1 <= order; ++l1) {
  //        for (int m1 = only_positive ? 0 : -l1; size_t(std::abs(m1)) <= l1; ++m1) {
  //          for (size_t l2 = 0; l2 <= order; ++l2) {
  //            for (int m2 = -int(l2); size_t(std::abs(m2)) <= l2; ++m2) {
  //              size_t row = pos(l1, m1);
  //              size_t col = pos(l2, only_positive ? std::abs(m2) : m2);
  //              RangeFieldType factor = !only_positive ? 1. : (m2 < 0 ? std::pow(-1., m2) : 1.);
  //              if (l1 == l2 + 1 && m1 == m2 + 1)
  //                By[row][col] += 0.5 * factor * std::complex<RangeFieldType>(0, 1) * B_lm(l2 + 1, m2 + 1);
  //              if (l1 == l2 - 1 && m1 == m2 + 1)
  //                By[row][col] += -0.5 * factor * std::complex<RangeFieldType>(0, 1) * B_lm(l2, -m2);
  //              if (l1 == l2 + 1 && m1 == m2 - 1)
  //                By[row][col] += 0.5 * factor * std::complex<RangeFieldType>(0, 1) * B_lm(l2 + 1, -m2 - 1);
  //              if (l1 == l2 - 1 && m1 == m2 - 1)
  //                By[row][col] += -0.5 * factor * std::complex<RangeFieldType>(0, 1) * B_lm(l2, m2);
  //            } // m2
  //          } // l2
  //        } // m1
  //      } // l1
  //      return By;
  //    } // ... create_By()

  static MatrixType create_Bz()
  {
    MatrixType Bz(dimRange, dimRange, 0);
    const auto& pos = helper<only_positive>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = only_positive ? 0. : -int(l1); size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          size_t row = pos(l1, m1);
          size_t col = pos(l2, m1); // m1 == m2, else matrix entry is 0
          if (l1 == l2 + 1)
            Bz[row][col] += A_lm(l2 + 1, m1);
          if (l1 == l2 - 1)
            Bz[row][col] += A_lm(l2, m1);
        } // l2
      } // m1
    } // l1
    return Bz;
  }

  template <bool positive, class anything = void>
  struct helper
  {
    // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
    // Each l has 2l+1 values of m, so (l, m) has position
    // (\sum_{k=0}^{l-1} (2k+1)) + (m+l) = l^2 + m + l
    static size_t pos(const size_t l, const int m)
    {
      return size_t(l * l + m + l);
    }
  };

  template <class anything>
  struct helper<true, anything>
  {
    // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
    // Each l has l+1 non-negative values of m, so (l, m) has position
    // (\sum_{k=0}^{l-1} (l+1)) + m = l(l+1)/2 + m
    static size_t pos(const size_t l, const int m)
    {
      return l * (l + 1) / 2 + m;
    }
  };
}; // class SphericalHarmonics<DomainFieldType, 3, ...>


template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_even = false>
class RealSphericalHarmonics
    : public BasisfunctionsInterface<DomainFieldType,
                                     3,
                                     RangeFieldType,
                                     only_even ? ((order + 1) * (order + 2)) / 2 : (order + 1) * (order + 1),
                                     1,
                                     fluxDim>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimFlux = fluxDim;
  static const size_t dimRange = only_even ? ((order + 1) * (order + 2)) / 2 : (order + 1) * (order + 1);

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, dimFlux> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  RealSphericalHarmonics(const QuadratureType& quadrature)
    : quadrature_(quadrature)
  {
  }

  virtual RangeType evaluate(const DomainType& v) const override
  {
    const auto v_spherical = XT::Common::CoordinateConverter<DomainFieldType>::to_spherical(v);
    return evaluate_in_spherical_coords(v_spherical);
  } // ... evaluate(...)

  RangeType evaluate_in_spherical_coords(const FieldVector<DomainFieldType, 2>& coords) const
  {
    const DomainFieldType theta = coords[0];
    const DomainFieldType phi = coords[1];
    RangeType ret(0);
    for (size_t ll = 0; ll <= order; ++ll)
      for (int mm = -int(ll); mm <= int(ll); ++mm)
        if (!only_even || !((mm + ll) % 2))
          ret[helper<only_even>::pos(ll, mm)] = evaluate_lm(theta, phi, int(ll), mm);
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override
  {
    RangeType ret(0);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 1;
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    return mass_matrix();
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimFlux> ret(MatrixType(dimRange, dimRange, 0));
    ret[0] = create_Bx();
    ret[1] = create_By();
    if (dimFlux == 3)
      ret[2] = create_Bz();
    return ret;
  } // ... mass_matrix_with_v()

  using BaseType::parallel_quadrature;

  // returns matrices with entries <v h_i h_j>_- and <v h_i h_j>_+
  virtual FieldVector<FieldVector<MatrixType, 2>, dimFlux> kinetic_flux_matrices() const
  {
    FieldVector<FieldVector<MatrixType, 2>, dimFlux> B_kinetic(
        FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    QuadratureType neg_quadrature;
    QuadratureType pos_quadrature;
    neg_quadrature.reserve(quadrature_.size());
    pos_quadrature.reserve(quadrature_.size());
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      neg_quadrature.clear();
      pos_quadrature.clear();
      for (const auto& quad_point : quadrature_) {
        const auto& v = quad_point.position();
        const auto& weight = quad_point.weight();
        if (XT::Common::FloatCmp::eq(v[dd], 0.)) {
          neg_quadrature.emplace_back(v, weight / 2.);
          pos_quadrature.emplace_back(v, weight / 2.);
        } else if (v[dd] > 0.)
          pos_quadrature.emplace_back(v, weight);
        else
          neg_quadrature.emplace_back(v, weight);
      }
      parallel_quadrature(neg_quadrature, B_kinetic[dd][0], dd);
      parallel_quadrature(pos_quadrature, B_kinetic[dd][1], dd);
    }
    return B_kinetic;
  } // ... kinetic_flux_matrices()

  virtual MatrixType reflection_matrix(const DomainType& n) const
  {
    MatrixType ret(dimRange, dimRange, 0);
    size_t direction;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (XT::Common::FloatCmp::ne(n[ii], 0.)) {
        direction = ii;
        if (XT::Common::FloatCmp::ne(std::abs(n[ii]), 1.))
          DUNE_THROW(NotImplemented, "Implemented only for +-e_i where e_i is the i-th canonical basis vector!");
      }
    }
    parallel_quadrature(quadrature_, ret, direction, true);
    ret.rightmultiply(mass_matrix_inverse());
    return ret;
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeType u_iso(0), alpha_iso(0);
    u_iso[0] = u[0];
    alpha_iso[0] = std::log(u[0] / std::sqrt(4. * M_PI)) * std::sqrt(4. * M_PI);
    return std::make_pair(u_iso, alpha_iso);
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      component_visualizer<DiscreteFunctionType, dimRange, 0>(u_n, filename_prefix, ii, std::sqrt(4 * M_PI));
    };
  }

  RangeFieldType calculate_psi_from_moments(const RangeType& val) const
  {
    return val[0] * std::sqrt(4 * M_PI);
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) { return XT::Common::to_string(val[0] * std::sqrt(4 * M_PI), 15); };
  } // ... stringifier()

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(u[0], u_bar[0]);
  }

private:
  static RangeFieldType A_lm(const size_t l, const int m)
  {
    return std::sqrt((l + m) * (l - m) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static RangeFieldType B_lm(const size_t l, const int m)
  {
    return std::sqrt((l + m) * (l + m - 1.) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static MatrixType create_Bx()
  {
    MatrixType Bx(dimRange, dimRange, 0.);
    const auto& pos = helper<only_even>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -int(l1); size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -int(l2); size_t(std::abs(m2)) <= l2; ++m2) {
            if (!only_even || (!((m1 + l1) % 2) && !((m2 + l2) % 2))) {
              if (l1 == l2 - 1 && m1 == m2 - 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += 0.5 * std::sqrt(1. + (m2 == 1)) * B_lm(l2, m2);
              if (l1 == l2 + 1 && m1 == m2 - 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += -0.5 * std::sqrt(1. + (m2 == 1)) * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == m2 + 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += -0.5 * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == m2 + 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += 0.5 * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == m2 + 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += 0.5 * (1. - (-m2 == 1)) * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == m2 + 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += -0.5 * (1. - (-m2 == 1)) * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == m2 - 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += -0.5 * B_lm(l2, m2);
              if (l1 == l2 + 1 && m1 == m2 - 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += 0.5 * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == 1 && m2 == 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += -1. / std::sqrt(2.) * B_lm(l2, 0);
              if (l1 == l2 + 1 && m1 == 1 && m2 == 0)
                Bx[pos(l1, m1)][pos(l2, m2)] += 1. / std::sqrt(2.) * B_lm(l2 + 1, 1);
            }
          } // m2
        } // l2
      } // m1
    } // l1
    return Bx;
  }

  static MatrixType create_By()
  {
    MatrixType By(dimRange, dimRange, 0.);
    const auto& pos = helper<only_even>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -int(l1); size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -int(l2); size_t(std::abs(m2)) <= l2; ++m2) {
            if (!only_even || (!((m1 + l1) % 2) && !((m2 + l2) % 2))) {
              if (l1 == l2 + 1 && m1 == -m2 + 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] += 0.5 * (1. - (m2 == 1)) * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == -m2 + 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] += -0.5 * (1. - (m2 == 1)) * B_lm(l2, m2);
              if (l1 == l2 - 1 && m1 == -m2 - 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] += -0.5 * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == -m2 - 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] += 0.5 * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == -m2 - 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] += 0.5 * std::sqrt(1. + (-m2 == 1)) * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == -m2 - 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] += -0.5 * std::sqrt(1. + (-m2 == 1)) * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == -m2 + 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] += 0.5 * B_lm(l2, m2);
              if (l1 == l2 + 1 && m1 == -m2 + 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] += -0.5 * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == -1 && m2 == 0)
                By[pos(l1, m1)][pos(l2, m2)] += -1. / std::sqrt(2.) * B_lm(l2, 0);
              if (l1 == l2 + 1 && m1 == -1 && m2 == 0)
                By[pos(l1, m1)][pos(l2, m2)] += 1. / std::sqrt(2.) * B_lm(l2 + 1, 1);
            }
          } // m2
        } // l2
      } // m1
    } // l1
    return By;
  } // ... create_By()

  static MatrixType create_Bz()
  {
    MatrixType Bz(dimRange, dimRange, 0);
    const auto& pos = helper<only_even>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -int(l1); size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -int(l2); size_t(std::abs(m2)) <= l2; ++m2) {
            if (!only_even || (!((m1 + l1) % 2) && !((m2 + l2) % 2))) {
              if (m1 == m2 && l1 == l2 + 1)
                Bz[pos(l1, m1)][pos(l2, m2)] += A_lm(l2 + 1, m2);
              if (m1 == m2 && l1 == l2 - 1)
                Bz[pos(l1, m1)][pos(l2, m2)] += A_lm(l2, m2);
            }
          } // m2
        } // l2
      } // m1
    } // l1
    return Bz;
  } // ... create_Bz()

  template <bool even, class anything = void>
  struct helper
  {
    // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
    // Each l has 2l+1 values of m, so (l, m) has position
    // (\sum_{k=0}^{l-1} (2k+1)) + (m+l) = l^2 + m + l
    static size_t pos(const size_t l, const int m)
    {
      return size_t(l * l + m + l);
    }
  };

  template <class anything>
  struct helper<true, anything>
  {
    // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
    // Each l has l+1 values of m (as only m s.t. m+l is even are considered), so (l, m) has position
    // (\sum_{k=0}^{l-1} (k+1)) + (m+l)/2 = l(l+1)/2 + (l+m)/2
    static size_t pos(const int l, const int m)
    {
      return size_t(l * (l + 1) / 2 + (m + l) / 2);
    }
  };

  // Notation from Garrett, Hauck, "A Comparison of Moment Closures for Linear Kinetic Transport Equations: The Line
  // Source Benchmark",
  // http://www.tandfonline.com/doi/full/10.1080/00411450.2014.910226?src=recsys&, Section 4.1
  RangeFieldType N_lm(const int l, const int m) const
  {
    assert(l >= 0 && m >= 0 && m <= l);
    return std::sqrt((2. * l + 1.) * XT::Common::factorial(l - m) / (XT::Common::factorial(l + m) * 4. * M_PI));
  }

  RangeFieldType evaluate_lm(const DomainFieldType theta, const DomainFieldType phi, const int l, const int m) const
  {
    const auto cos_theta = std::cos(theta);
    assert(l >= 0 && std::abs(m) <= l);
    if (m < 0)
      return std::sqrt(2) * N_lm(l, -m) * boost::math::legendre_p(l, -m, cos_theta) * std::sin(-m * phi);
    else if (m == 0)
      return N_lm(l, 0) * boost::math::legendre_p(l, 0, cos_theta);
    else
      return std::sqrt(2) * N_lm(l, m) * boost::math::legendre_p(l, m, cos_theta) * std::cos(m * phi);
  }

  const QuadratureType& quadrature_;
}; // class RealSphericalHarmonics<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_SPHERICALHARMONICS_HH
