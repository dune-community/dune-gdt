// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH

#include <memory>
#include <vector>
#include <string>

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/triangulation.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/lebedevquadrature.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


// see https://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion
template <class FieldType, int rows>
Dune::FieldMatrix<FieldType, rows, rows> tridiagonal_matrix_inverse(const FieldMatrix<FieldType, rows, rows>& matrix)
{
  typedef Dune::FieldMatrix<FieldType, rows, rows> MatrixType;
  size_t cols = rows;
#ifndef NDEBUG
  for (size_t rr = 0; rr < rows; ++rr)
    for (size_t cc = 0; cc < cols; ++cc)
      if ((cc > rr + 1 || cc + 1 < rr) && XT::Common::FloatCmp::ne(matrix[rr][cc], 0.))
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "Matrix has to be tridiagonal!");
#endif // NDEBUG
  MatrixType ret(0);
  Dune::FieldVector<FieldType, rows + 1> a(0), b(0), c(0), theta(0);
  Dune::FieldVector<FieldType, rows + 2> phi(0);
  for (size_t ii = 1; ii < rows + 1; ++ii) {
    a[ii] = matrix[ii - 1][ii - 1];
    if (ii < rows) {
      b[ii] = matrix[ii - 1][ii];
      c[ii] = matrix[ii][ii - 1];
    }
  }
  theta[0] = 1;
  theta[1] = a[1];
  for (size_t ii = 2; ii < rows + 1; ++ii)
    theta[ii] = a[ii] * theta[ii - 1] - b[ii - 1] * c[ii - 1] * theta[ii - 2];
  phi[rows + 1] = 1;
  phi[rows] = a[rows];
  for (size_t ii = rows - 1; ii > 0; --ii)
    phi[ii] = a[ii] * phi[ii + 1] - b[ii] * c[ii] * phi[ii + 2];
  for (size_t ii = 1; ii < rows + 1; ++ii) {
    for (size_t jj = 1; jj < cols + 1; ++jj) {
      if (ii == jj)
        ret[ii - 1][jj - 1] = theta[ii - 1] * phi[jj + 1] / theta[rows];
      else if (ii < jj) {
        ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[ii - 1] * phi[jj + 1] / theta[rows];
        for (size_t kk = ii; kk < jj; ++kk)
          ret[ii - 1][jj - 1] *= b[kk];
      } else if (ii > jj) {
        ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[jj - 1] * phi[ii + 1] / theta[rows];
        for (size_t kk = jj; kk < ii; ++kk)
          ret[ii - 1][jj - 1] *= c[kk];
      }
    } // jj
  } // ii
  return ret;
} // ... tridiagonal_matrix_inverse(...)


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain>
class BasisfunctionsInterface
{
public:
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;
  typedef typename XT::Functions::RangeTypeSelector<RangeFieldType, dimRange, dimRangeCols>::type RangeType;

  virtual ~BasisfunctionsInterface(){};

  virtual RangeType evaluate(const DomainType& v) const = 0;

  virtual RangeType integrated() const = 0;

  virtual MatrixType mass_matrix() const = 0;

  virtual MatrixType mass_matrix_inverse() const = 0;

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const = 0;
};


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain>
class HatFunctions
{
  //  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain,
          size_t order = 1>
class PiecewiseMonomials
{
  //  static_assert(false, "Not implemented for this combination of dimension and order!");
};

template <class DomainFieldType, class RangeFieldType, size_t order, size_t dimRangeCols = 1>
class LegendrePolynomials : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, order + 1, dimRangeCols>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = order + 1;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;

  virtual RangeType evaluate(const DomainType& v) const override
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      for (size_t kk = 0; kk <= ii; ++kk)
        ret[ii] += std::pow(-v, kk) * XT::Common::binomial_coefficient(ii, kk)
                   * XT::Common::binomial_coefficient((ii + kk - 1) / 2., ii);
      ret[ii] *= std::pow(-1, ii) * (1 << ii); // (-2)^ii
    }
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override
  {
    RangeType ret(0);
    ret[0] = 2;
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 2. / (2. * rr + 1.);
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    MatrixType Minv(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      Minv[rr][rr] = (2. * rr + 1.) / 2.;
    return Minv;
  }

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    MatrixType B(0);
    for (size_t rr = 0; rr < dimRange; ++rr) {
      for (size_t cc = 0; cc < dimRange; ++cc) {
        if (cc == rr - 1)
          B[rr][cc] = 2 * rr / (4. * rr * rr - 1.);
        else if (cc == rr + 1)
          B[rr][cc] = (2 * rr + 2.) / ((2. * rr + 1.) * (2. * rr + 3));
      }
    }
    return B;
  }

  MatrixType S() const
  {
    MatrixType S(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      S[rr][rr] = -2. * rr * (rr + 1.) / (2 * rr + 1);
    return S;
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeType u_iso(0), alpha_iso(0);
    u_iso[0] = u[0];
    alpha_iso[0] = std::log(u[0] / 2.);
    return std::make_pair(u_iso, alpha_iso);
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(u[0], u_bar[0]);
  }
}; // class LegendrePolynomials<DomainFieldType, 1, ...>

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_positive = false>
class SphericalHarmonics
    : public BasisfunctionsInterface<DomainFieldType,
                                     3,
                                     RangeFieldType,
                                     only_positive ? ((order + 1) * (order + 2)) / 2 : (order + 1) * (order + 1),
                                     1,
                                     fluxDim>
{
  static const size_t dimDomain = 3;
  static const size_t dimRange = only_positive ? ((order + 1) * (order + 2)) / 2 : (order + 1) * (order + 1);
  static const size_t dimFlux = fluxDim;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, dimFlux> BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;

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
    for (int ll = 0; ll <= order; ++ll)
      for (int mm = only_positive ? 0 : -ll; mm <= ll; ++mm)
        ret[helper<only_positive>::pos(ll, mm)] = boost::math::spherical_harmonic(ll, mm, theta, phi);
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
    MatrixType M(0);
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
    FieldVector<MatrixType, dimFlux> ret(MatrixType(0));
    ret[0] = create_Bx();
    ret[1] = create_By();
    if (dimFlux == 3)
      ret[2] = create_Bz();
    return ret;
  } // ... mass_matrix_with_v()

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
    MatrixType Bx(0);
    const auto& pos = helper<only_positive>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = only_positive ? 0 : -l1; std::abs(m1) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = only_positive ? 0 : -l2; std::abs(m2) <= l2; ++m2) {
            if (l1 == l2 + 1 && m1 == m2 + 1)
              Bx[pos(l1, m1)][pos(l2, m2)] = -0.5 * B_lm(l2 + 1, m2 + 1);
            if (l1 == l2 - 1 && m1 == m2 + 1)
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * B_lm(l2, -m2);
            if (l1 == l2 + 1 && m1 == m2 - 1)
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * B_lm(l2 + 1, -m2 - 1);
            if (l1 == l2 - 1 && m1 == m2 - 1)
              Bx[pos(l1, m1)][pos(l2, m2)] = -0.5 * B_lm(l2, m2);
          } // m2
        } // l2
      } // m1
    } // l1
    return Bx;
  } // ... create_Bx()

  static MatrixType create_By()
  {
    MatrixType By(0);
    const auto& pos = helper<only_positive>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = only_positive ? 0 : -l1; std::abs(m1) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = only_positive ? 0 : -l2; std::abs(m2) <= l2; ++m2) {
            if (l1 == l2 + 1 && m1 == m2 + 1)
              By[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::complex<RangeFieldType>(0, 1) * B_lm(l2 + 1, m2 + 1);
            if (l1 == l2 - 1 && m1 == m2 + 1)
              By[pos(l1, m1)][pos(l2, m2)] = -0.5 * std::complex<RangeFieldType>(0, 1) * B_lm(l2, -m2);
            if (l1 == l2 + 1 && m1 == m2 - 1)
              By[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::complex<RangeFieldType>(0, 1) * B_lm(l2 + 1, -m2 - 1);
            if (l1 == l2 - 1 && m1 == m2 - 1)
              By[pos(l1, m1)][pos(l2, m2)] = -0.5 * std::complex<RangeFieldType>(0, 1) * B_lm(l2, m2);
          } // m2
        } // l2
      } // m1
    } // l1
    return By;
  } // ... create_By()

  static MatrixType create_Bz()
  {
    MatrixType Bz(0);
    const auto& pos = helper<only_positive>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = only_positive ? 0 : -l1; std::abs(m1) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = only_positive ? 0 : -l2; std::abs(m2) <= l2; ++m2) {
            if (m1 == m2 && l1 == l2 + 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = A_lm(l2 + 1, m2);
            if (m1 == m2 && l1 == l2 - 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = A_lm(l2, m2);
          } // m2
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
    static size_t pos(const int l, const int m)
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

template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

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
          ret[helper<only_even>::pos(ll, mm)] = evaluate_lm(theta, phi, ll, mm);
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
    MatrixType M(0);
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
    FieldVector<MatrixType, dimFlux> ret(MatrixType(0));
    ret[0] = create_Bx();
    ret[1] = create_By();
    if (dimFlux == 3)
      ret[2] = create_Bz();
    return ret;
  } // ... mass_matrix_with_v()

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeType u_iso(0), alpha_iso(0);
    u_iso[0] = u[0];
    alpha_iso[0] = std::log(u[0] / std::sqrt(4. * M_PI)) * std::sqrt(4. * M_PI);
    return std::make_pair(u_iso, alpha_iso);
  }

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
    MatrixType Bx(0);
    const auto& pos = helper<only_even>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -l1; size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -l2; size_t(std::abs(m2)) <= l2; ++m2) {
            if (!only_even || (!((m1 + l1) % 2) && !((m2 + l2) % 2))) {
              if (l1 == l2 - 1 && m1 == m2 - 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt(1. + (m2 == 1)) * B_lm(l2, m2);
              if (l1 == l2 + 1 && m1 == m2 - 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = -0.5 * std::sqrt(1. + (m2 == 1)) * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == m2 + 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = -0.5 * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == m2 + 1 && m2 > 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == m2 + 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * (1. - (-m2 == 1)) * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == m2 + 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = -0.5 * (1. - (-m2 == 1)) * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == m2 - 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = -0.5 * B_lm(l2, m2);
              if (l1 == l2 + 1 && m1 == m2 - 1 && m2 < 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == 1 && m2 == 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = -1. / std::sqrt(2.) * B_lm(l2, 0);
              if (l1 == l2 + 1 && m1 == 1 && m2 == 0)
                Bx[pos(l1, m1)][pos(l2, m2)] = 1. / std::sqrt(2.) * B_lm(l2 + 1, 1);
            }
          } // m2
        } // l2
      } // m1
    } // l1
    return Bx;
  }

  static MatrixType create_By()
  {
    MatrixType By(0);
    const auto& pos = helper<only_even>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -l1; size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -l2; size_t(std::abs(m2)) <= l2; ++m2) {
            if (!only_even || (!((m1 + l1) % 2) && !((m2 + l2) % 2))) {
              if (l1 == l2 + 1 && m1 == -m2 + 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] = 0.5 * (1. - (m2 == 1)) * B_lm(l2 + 2, -m2 + 1);
              if (l1 == l2 - 1 && m1 == -m2 + 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] = -0.5 * (1. - (m2 == 1)) * B_lm(l2, m2);
              if (l1 == l2 - 1 && m1 == -m2 - 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] = -0.5 * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == -m2 - 1 && m2 > 0)
                By[pos(l1, m1)][pos(l2, m2)] = 0.5 * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == -m2 - 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt(1. + (-m2 == 1)) * B_lm(l2, -m2);
              if (l1 == l2 + 1 && m1 == -m2 - 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] = -0.5 * std::sqrt(1. + (-m2 == 1)) * B_lm(l2 + 1, m2 + 1);
              if (l1 == l2 - 1 && m1 == -m2 + 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] = 0.5 * B_lm(l2, m2);
              if (l1 == l2 + 1 && m1 == -m2 + 1 && m2 < 0)
                By[pos(l1, m1)][pos(l2, m2)] = -0.5 * B_lm(l2 + 1, -m2 + 1);
              if (l1 == l2 - 1 && m1 == -1 && m2 == 0)
                By[pos(l1, m1)][pos(l2, m2)] = -1. / std::sqrt(2.) * B_lm(l2, 0);
              if (l1 == l2 + 1 && m1 == -1 && m2 == 0)
                By[pos(l1, m1)][pos(l2, m2)] = 1. / std::sqrt(2.) * B_lm(l2 + 1, 1);
            }
          } // m2
        } // l2
      } // m1
    } // l1
    return By;
  } // ... create_By()

  static MatrixType create_Bz()
  {
    MatrixType Bz(0);
    const auto& pos = helper<only_even>::pos;
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -l1; size_t(std::abs(m1)) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -l2; size_t(std::abs(m2)) <= l2; ++m2) {
            if (!only_even || (!((m1 + l1) % 2) && !((m2 + l2) % 2))) {
              if (m1 == m2 && l1 == l2 + 1)
                Bz[pos(l1, m1)][pos(l2, m2)] = A_lm(l2 + 1, m2);
              if (m1 == m2 && l1 == l2 - 1)
                Bz[pos(l1, m1)][pos(l2, m2)] = A_lm(l2, m2);
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
}; // class RealSphericalHarmonics<DomainFieldType, 3, ...>

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols>
class HatFunctions<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  typedef RangeType TriangulationType;

  HatFunctions(const TriangulationType triangulation = create_triangulation(),
               const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    RangeType ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = -1. + 2. * ii / (dimRange - 1.);
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      if (ii < dimRange - 1 && XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1]))
        ret[ii] = (v - triangulation_[ii + 1]) / (triangulation_[ii] - triangulation_[ii + 1]);
      if (ii > 0 && XT::Common::FloatCmp::ge(v[0], triangulation_[ii - 1])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii]))
        ret[ii] = (v - triangulation_[ii - 1]) / (triangulation_[ii] - triangulation_[ii - 1]);
    }
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override
  {
    RangeType ret(0);
    ret[0] = triangulation_[1] - triangulation_[0];
    for (size_t ii = 1; ii < dimRange - 1; ++ii)
      ret[ii] = triangulation_[ii + 1] - triangulation_[ii - 1];
    ret[dimRange - 1] = triangulation_[dimRange - 1] - triangulation_[dimRange - 2];
    ret *= 0.5;
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override
  {
    MatrixType ret(0);
    ret[0][0] = (triangulation_[1] - triangulation_[0]) / 3.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < dimRange - 1)
        ret[rr][rr] = (triangulation_[rr + 1] - triangulation_[rr - 1]) / 3.;
      if (rr > 0)
        ret[rr][rr - 1] = (triangulation_[rr] - triangulation_[rr - 1]) / 6.;
      if (rr < dimRange - 1)
        ret[rr][rr + 1] = (triangulation_[rr + 1] - triangulation_[rr]) / 6.;
    }
    ret[dimRange - 1][dimRange - 1] = (triangulation_[dimRange - 1] - triangulation_[dimRange - 2]) / 3.;
    return ret;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    return tridiagonal_matrix_inverse(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, 1> mass_matrix_with_v() const override
  {
    MatrixType ret(0);
    ret[0][0] = (triangulation_[1] * triangulation_[1] + 2 * triangulation_[1] * triangulation_[0]
                 - 3 * triangulation_[0] * triangulation_[0])
                / 12.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < dimRange - 1)
        ret[rr][rr] = (triangulation_[rr + 1] * triangulation_[rr + 1] + 2 * triangulation_[rr + 1] * triangulation_[rr]
                       - 2 * triangulation_[rr] * triangulation_[rr - 1]
                       - triangulation_[rr - 1] * triangulation_[rr - 1])
                      / 12.;
      if (rr > 0)
        ret[rr][rr - 1] =
            (triangulation_[rr] * triangulation_[rr] - triangulation_[rr - 1] * triangulation_[rr - 1]) / 12.;
      if (rr < dimRange - 1)
        ret[rr][rr + 1] =
            (triangulation_[rr + 1] * triangulation_[rr + 1] - triangulation_[rr] * triangulation_[rr]) / 12.;
    }
    ret[dimRange - 1][dimRange - 1] = (3 * triangulation_[dimRange - 1] * triangulation_[dimRange - 1]
                                       - 2 * triangulation_[dimRange - 1] * triangulation_[dimRange - 2]
                                       - triangulation_[dimRange - 2] * triangulation_[dimRange - 2])
                                      / 12.;
    return ret;
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      psi_iso += u[ii];
    psi_iso /= 2.;
    RangeType alpha_iso(std::log(psi_iso)), u_iso;
    u_iso = integrated();
    u_iso *= psi_iso / 2.;
    return std::make_pair(u_iso, alpha_iso);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(std::accumulate(u.begin(), u.end(), RangeFieldType(0)),
                        std::accumulate(u_bar.begin(), u_bar.end(), RangeFieldType(0)));
  }

private:
  const TriangulationType triangulation_;
}; // class HatFunctions<DomainFieldType, 1, ...>


// After each refinement step:
// num_vertices_new = num_vertices_old + num_intersections_old
// num_intersections_new = 2*num_intersections_old + 3*num_faces_old
// num_faces_new = 4*num_faces_old
// Initially, there are 6 vertices, 12 intersections and 8 faces.
template <size_t refinements>
struct OctaederStatistics
{
  static constexpr size_t num_faces()
  {
    return 8 * (1 << 2 * refinements);
  }

  static constexpr size_t num_intersections()
  {
    return 2 * OctaederStatistics<refinements - 1>::num_intersections()
           + 2 * OctaederStatistics<refinements - 1>::num_faces();
  }

  static constexpr size_t num_vertices()
  {
    return OctaederStatistics<refinements - 1>::num_vertices()
           + OctaederStatistics<refinements - 1>::num_intersections();
  }
};

template <>
struct OctaederStatistics<0>
{
  static constexpr size_t num_faces()
  {
    return 8;
  }

  static constexpr size_t num_intersections()
  {
    return 12;
  }

  static constexpr size_t num_vertices()
  {
    return 6;
  }
};

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t fluxDim>
class HatFunctions<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static const size_t dimFlux = fluxDim;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, dimFlux> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef SphericalTriangulation<DomainFieldType> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;

  HatFunctions(const TriangulationType& triangulation, const QuadratureType& quadrature)
    : triangulation_(triangulation)
    , quadrature_(quadrature)
  {
    assert(triangulation_.vertices().size() == dimRange);
  }

  HatFunctions(const size_t refinements = 0,
               const size_t quadrature_refinements = 4,
               std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                   {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements)
    , quadrature_(triangulation_.quadrature_rule(quadrature_refinements))
  {
    assert(triangulation_.vertices().size() == dimRange);
  }

  virtual RangeType evaluate(const DomainType& v) const override
  {
    RangeType ret(0);
    bool success = false;
    // walk over faces
    for (const auto& face : triangulation_.faces()) {
      const auto& vertices = face->vertices();
      DomainType barycentric_coords(0);
      success = calculate_barycentric_coordinates(v, vertices, barycentric_coords);
      if (success) {
        for (size_t ii = 0; ii < 3; ++ii) {
          ret[vertices[ii]->index()] = barycentric_coords[ii];
        }
        std::cout << "vertices: " << XT::Common::to_string(face->vertices()[0]->position()) << ", "
                  << XT::Common::to_string(face->vertices()[1]->position()) << ", "
                  << XT::Common::to_string(face->vertices()[2]->position())
                  << ", coords: " << XT::Common::to_string(ret) << ", v: " << XT::Common::to_string(v) << std::endl;
        break;
      }
    } // faces
    assert(success);
    return ret;
  } // ... evaluate(...)

  // avoid recalculation of integral by using a static local variable that is initialized on first call
  virtual RangeType integrated() const override
  {
    static const RangeType ret = integrated_initializer();
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType A(0);
    for (const auto& quad_point : quadrature_) {
      const auto basis_evaluated = evaluate(quad_point.position());
      for (size_t nn = 0; nn < dimRange; ++nn) {
        for (size_t mm = 0; mm < dimRange; ++mm) {
          A[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * quad_point.weight();
        } // mm
      } // nn
    } // quadrature
    return A;
  } // ... create_flux_config(...)

  virtual MatrixType mass_matrix_inverse() const override
  {
    auto ret = mass_matrix();
    ret.invert();
    return ret;
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(0));
    for (const auto& quad_point : quadrature_) {
      const auto& v = quad_point.position();
      const auto basis_evaluated = evaluate(v);
      const auto& weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        for (size_t mm = 0; mm < dimRange; ++mm)
          for (size_t dd = 0; dd < dimFlux; ++dd)
            B[dd][nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * v[dd] * weight;
    } // quadrature
    return B;
  } // ... mass_matrix_with_v()

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      psi_iso += u[ii];
    psi_iso /= 4. * M_PI;
    RangeType alpha_iso(std::log(psi_iso));
    auto u_iso = integrated();
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(std::accumulate(u.begin(), u.end(), RangeFieldType(0)),
                        std::accumulate(u_bar.begin(), u_bar.end(), RangeFieldType(0)));
  }

  const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

protected:
  RangeType integrated_initializer() const
  {
    RangeType ret(0);
    for (const auto& quad_point : quadrature_) {
      const auto basis_evaluated = evaluate(quad_point.position());
      basis_evaluated *= quad_point.weight();
      ret += basis_evaluated;
    } // quadrature
    return ret;
  }

  template <class VertexVectorType>
  bool calculate_barycentric_coordinates(const DomainType& v, const VertexVectorType& vertices, DomainType& ret) const
  {
    Dune::FieldMatrix<RangeFieldType, 3, 3> gradients(0);
    for (size_t ii = 0; ii < 3; ++ii) {
      // copy vertices to gradients
      gradients[ii] = vertices[ii]->position();
      const auto scalar_prod = v * gradients[ii];
      // if v is not on the same octant of the sphere as the vertices, return false
      // assumes the triangulation is fine enough that vertices[ii]*vertices[jj] >= 0 for all triangles
      if (XT::Common::FloatCmp::lt(scalar_prod, 0.))
        return false;
      auto v_scaled = v;
      v_scaled *= scalar_prod;
      gradients[ii] -= v_scaled;
      // scale with factor
      auto denominator = std::sqrt(1. - std::pow(scalar_prod, 2));
      gradients[ii] *= XT::Common::FloatCmp::eq(denominator, 0.) ? 0. : std::acos(scalar_prod) / denominator;
    } // ii
    // Calculate barycentric coordinates for 0 w.r.t to the points g_i = gradients[i]
    // For that purpose, solve the overdetermined system  A (h0 h1)^T = b
    // for the matrix A = (g_0-g_2 g_1-g_2) and the right-hand side b = -g_2.
    // The solution is (A^T A)^{-1} A^T b.
    // The third coordinate is calculated from the condition h0+h1+h2=1.
    Dune::FieldMatrix<RangeFieldType, 3, 2> A;
    Dune::FieldMatrix<RangeFieldType, 2, 3> AT;
    Dune::FieldVector<RangeFieldType, 2> solution;
    AT[0] = gradients[0];
    AT[1] = gradients[1];
    AT[0] -= gradients[2];
    AT[1] -= gradients[2];
    for (size_t ii = 0; ii < 3; ++ii)
      for (size_t jj = 0; jj < 2; ++jj)
        A[ii][jj] = AT[jj][ii];
    Dune::FieldMatrix<RangeFieldType, 2, 2> AT_A = AT.rightmultiplyany(A);
    gradients[2] *= -1;
    FieldVector<RangeFieldType, 2> AT_b;
    AT.mv(gradients[2], AT_b);
    AT_A.solve(solution, AT_b);
    ret[0] = solution[0];
    ret[1] = solution[1];
    ret[2] = 1. - ret[0] - ret[1];
    if (XT::Common::FloatCmp::lt(ret[0], 0.) || XT::Common::FloatCmp::lt(ret[1], 0.))
      return false;
    if (XT::Common::FloatCmp::lt(ret[2], 0.))
      return false;
    return true;
  } // bool calculate_barycentric_coordinates(...)

  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
}; // class HatFunctions<DomainFieldType, 3, ...>


template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t dimFlux>
class PiecewiseMonomials<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static_assert(!(dimRange % 2), "dimRange has to be even!");

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef FieldVector<DomainFieldType, dimRange + 1> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;

  PiecewiseMonomials(const TriangulationType& triangulation = create_triangulation(),
                     const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    TriangulationType ret;
    for (size_t ii = 0; ii < dimRange / 2 + 1; ++ii)
      ret[ii] = -1. + 4. * ii / dimRange;
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      if (XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1])) {
        ret[2 * ii] = 1;
        ret[2 * ii + 1] = v[0];
      }
    }
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      ret[2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      ret[2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
    }
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      M[2 * ii][2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      M[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      M[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      M[2 * ii + 1][2 * ii] = M[2 * ii][2 * ii + 1];
    }
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    return tridiagonal_matrix_inverse(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    MatrixType B(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      B[2 * ii][2 * ii] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      B[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 4) - std::pow(triangulation_[ii], 4)) / 4.;
      B[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      B[2 * ii + 1][2 * ii] = B[2 * ii][2 * ii + 1];
    }
    return FieldVector<MatrixType, dimDomain>(B);
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeType alpha_iso(0);
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ii += 2) {
      psi_iso += u[ii];
      alpha_iso[ii] = 1;
    }
    psi_iso /= 2.;
    alpha_iso *= std::log(psi_iso);
    RangeType u_iso = integrated();
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    RangeFieldType u_sum;
    auto u_bar_sum = u_sum;
    for (size_t ii = 0; ii < u.size(); ii += 4) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
  }

private:
  const TriangulationType triangulation_;
}; // class PiecewiseMonomials<DomainFieldType, 1, ...>


template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t dimFlux>
class PiecewiseMonomials<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, dimFlux>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, dimFlux> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef SphericalTriangulation<DomainFieldType> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;

  PiecewiseMonomials(const TriangulationType& triangulation, const QuadratureType& quadrature)
    : triangulation_(triangulation)
    , quadrature_(quadrature)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  PiecewiseMonomials(const size_t refinements = 0,
                     const size_t quadrature_refinements = 4,
                     std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                         {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements)
    , quadrature_(triangulation_.quadrature_rule(quadrature_refinements))
  {
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    RangeType ret(0);
    FieldMatrix<RangeFieldType, 3, 3> vertices_matrix;
    FieldMatrix<RangeFieldType, 3, 3> determinant_matrix;
    for (const auto& face : triangulation_.faces()) {
      // vertices are ordered counterclockwise, so if the points is inside the spherical triangle,
      // the coordinate system formed by two adjacent vertices and v is always right-handed, i.e.
      // the triple product is positive
      const auto& vertices = face->vertices();
      for (size_t ii = 0; ii < 3; ++ii)
        vertices_matrix[ii] = vertices[ii]->position();
      bool v_in_this_facet = true;
      // the triple products that need to be positive are the determinants of the matrices (v1, v2, v), (v2, v3, v),
      // (v3, v1, v), where vi is the ith vertex. Swapping two columns changes the sign of det, the matrices used
      // below all have an even number of column swaps
      for (size_t ii = 0; ii < 3; ++ii) {
        determinant_matrix = vertices_matrix;
        determinant_matrix[ii] = v;
        if (XT::Common::FloatCmp::lt(determinant_matrix.determinant(), 0.)) {
          v_in_this_facet = false;
          break;
        }
      }
      if (v_in_this_facet) {
        const auto face_index = face->index();
        ret[4 * face_index] = 1;
        for (size_t ii = 1; ii < 4; ++ii)
          ret[4 * face_index + ii] = v[ii - 1];
        break;
      }
    } // faces
    return ret;
  } // ... evaluate(...)

  // returns <b>, where b is the basis functions vector
  virtual RangeType integrated() const override
  {
    static const RangeType ret = integrated_initializer();
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(0);
    for (const auto& quad_point : quadrature_) {
      const auto basis_evaluated = evaluate(quad_point.position());
      for (size_t nn = 0; nn < dimRange; ++nn)
        for (size_t mm = 0; mm < dimRange; ++mm)
          M[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * quad_point.weight();
    } // quadrature
    return M;
  } // ... mass_matrix()

  virtual MatrixType mass_matrix_inverse() const override
  {
    auto ret = mass_matrix();
    ret.invert();
    return ret;
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(0));
    for (const auto& quad_point : quadrature_) {
      const auto v = quad_point.position();
      const auto basis_evaluated = evaluate(v);
      const auto weight = quad_point.weight();
      for (size_t dd = 0; dd < dimFlux; ++dd)
        for (size_t nn = 0; nn < dimRange; ++nn)
          for (size_t mm = 0; mm < dimRange; ++mm)
            B[dd][nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * v[dd] * weight;
    } // quadrature
    return B;
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    RangeType alpha_iso(0);
    for (size_t ii = 0; ii < dimRange; ii += 4) {
      psi_iso += u[ii];
      alpha_iso[ii] = 1.;
    }
    psi_iso /= 4. * M_PI;
    alpha_iso *= std::log(psi_iso);
    auto u_iso = integrated();
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    RangeFieldType u_sum;
    auto u_bar_sum = u_sum;
    for (size_t ii = 0; ii < u.size(); ii += 4) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

private:
  RangeType integrated_initializer() const
  {
    RangeType ret(0);
    for (const auto& quad_point : quadrature_) {
      auto basis_evaluated = evaluate(quad_point.position());
      basis_evaluated *= quad_point.weight();
      ret += basis_evaluated;
    } // quadrature
    return ret;
  }

  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
}; // class PiecewiseMonomials<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH
