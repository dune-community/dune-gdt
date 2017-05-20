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

#include <boost/geometry.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/la/container.hh>

#include "../default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

// see https://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion
template <class MatrixType>
MatrixType tridiagonal_matrix_inverse(const MatrixType& matrix)
{
#ifndef NDEBUG
  for (size_t rr = 0; rr < dimRange; ++rr)
    for (size_t cc = 0; cc < dimRange; ++cc)
      if ((cc > rr + 1 || cc + 1 < rr) && XT::Common::FloatCmp::ne(matrix[rr][cc], 0.))
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "Matrix has to be tridiagonal!");
#endif // NDEBUG
  MatrixType ret(0);
  Dune::FieldVector<RangeFieldType, dimRange + 1> a(0), b(0), c(0), theta(0);
  Dune::FieldVector<RangeFieldType, dimRange + 2> phi(0);
  for (size_t ii = 1; ii < dimRange + 1; ++ii) {
    a[ii] = matrix[ii - 1][ii - 1];
    if (ii < dimRange) {
      b[ii] = matrix[ii - 1][ii];
      c[ii] = matrix[ii][ii - 1];
    }
  }
  theta[0] = 1;
  theta[1] = a[1];
  for (size_t ii = 2; ii < dimRange + 1; ++ii)
    theta[ii] = a[ii] * theta[ii - 1] - b[ii - 1] * c[ii - 1] * theta[ii - 2];
  phi[dimRange + 1] = 1;
  phi[dimRange] = a[dimRange];
  for (size_t ii = dimRange - 1; ii > 0; --ii)
    phi[ii] = a[ii] * phi[ii + 1] - b[ii] * c[ii] * phi[ii + 2];
  for (size_t ii = 1; ii < dimRange + 1; ++ii) {
    for (size_t jj = 1; jj < dimRange + 1; ++jj) {
      if (ii == jj)
        ret[ii - 1][jj - 1] = theta[ii - 1] * phi[jj + 1] / theta[dimRange];
      else if (ii < jj) {
        ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[ii - 1] * phi[jj + 1] / theta[dimRange];
        for (size_t kk = ii; kk < jj; ++kk)
          ret[ii - 1][jj - 1] *= b[kk];
      } else if (ii > jj) {
        ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[jj - 1] * phi[ii + 1] / theta[dimRange];
        for (size_t kk = jj; kk < ii; ++kk)
          ret[ii - 1][jj - 1] *= c[kk];
      }
    } // jj
  } // ii
  return ret;
} // ... tridiagonal_matrix_inverse(...)


template <class DomainFieldType, size_t dimDomain, class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1>
class BasisfunctionsInterface
{
public:
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename XT::Functions::internal::RangeTypeSelector<RangeFieldType, dimRange, dimRangeCols>::type RangeType;

  virtual RangeType evaluate(const DomainType& v) const = 0;

  virtual RangeType integrated() const = 0;

  virtual MatrixType mass_matrix() const = 0;

  virtual MatrixType mass_matrix_inverse() const = 0;

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const = 0;
};


template <class DomainFieldType, size_t dimDomain, class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1>
class LegendrePolynomials
{
  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType, size_t dimDomain, class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1>
class RealSphericalHarmonics
{
  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType, size_t dimDomain, class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1>
class HatFunctions
{
  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t order = 1>
class PiecewiseMonomials
{
  static_assert(false, "Not implemented for this combination of dimension and order!");
};

template <class DomainFieldType, class RangeFieldType, size_t dimRange, size_t dimRangeCols>
class LegendrePolynomials<DomainFieldType, 1, RangeFieldType, dimRange, dimRangeCols>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, dimRange, dimRangeCols>
{
  static const size_t dimDomain = 1;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using BaseType::DomainType;
  using BaseType::RangeType;

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

  virtual FieldVector<MatrixType, dimRange> mass_matrix_with_v() const override
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
}; // class LegendrePolynomials<DomainFieldType, 1, ...>

// converts from (x, y, z) to (theta, phi) on the unit sphere s.t.
// (x, y, z) = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta))
template <class DomainFieldType>
struct CoordinateConverter
{
  typedef typename FieldVector<DomainFieldType, 3> CartesianCoordType;
  typedef typename FieldVector<DomainFieldType, 2> SphericalCoordType;
  typedef typename boost::geometry::point<DomainFieldType, dimDomain, boost::geometry::cs::cartesian>
      BoostCartesianCoordType;
  typedef typename boost::geometry::point<DomainFieldType, 2, boost::geometry::cs::spherical<boost::geometry::radian>>
      BoostSphericalCoordType;

  static SphericalCoordType to_spherical(const CartesianCoordType& x)
  {
    BoostCartesianCoordType x_boost(x[0], x[1], x[2]);
    BoostSphericalCoordType x_spherical_boost;
    boost::geometry::transform(x_boost, x_spherical_boost);
    return SphericalCoordType{x_spherical_boost[0], x_spherical_boost[1]};
  }

  static CartesianCoordType to_cartesian(const SphericalCoordType& x_spherical)
  {
    BoostSphericalCoordType x_spherical_boost(x_spherical[0], x_spherical[1]);
    BoostCartesianCoordType x_boost;
    boost::geometry::transform(x_spherical_boost, x_boost);
    return CartesianCoordType{x_boost[0], x_boost[1], x_boost[2]};
  }
};


template <class DomainFieldType, class RangeFieldType>
struct RealSphericalHarmonicsBase<DomainFieldType, 3, RangeFieldType>
{
  static RangeFieldType evaluate_associated_legendre_polynomial(const FieldType& mu, const int l, int m) const
  {
    return boost::math::legendre_p(l, m, mu);
  }

  // Notation from Garrett, Hauck, "A Comparison of Moment Closures for Linear Kinetic Transport Equations: The Line
  // Source Benchmark",
  // http://www.tandfonline.com/doi/full/10.1080/00411450.2014.910226?src=recsys&, Section 4.1
  RangeFieldType N_lm(const int l, const int m) const
  {
    assert(l >= 0 && m >= 0 && m <= l);
    return std::sqrt((2. * l + 1.) * XT::Common::factorial(l - m) / (XT::Common::factorial(l + m) * 4. * M_PI));
  }

  template <class FieldType>
  DomainFieldType evaluate_real_spherical_harmonics(const DomainFieldType theta,
                                                    const DomainFieldType phi,
                                                    const int l,
                                                    const int m) const
  {
    const auto cos_theta = std::cos(theta);
    assert(l >= 0 && std::abs(m) <= l);
    if (m < 0)
      return std::sqrt(2) * N_lm(l, -m) * evaluate_associated_legendre_polynomial(cos_theta, l, -m)
             * std::sin(-m * phi);
    else if (m == 0)
      return N_lm(l, 0) * evaluate_associated_legendre_polynomial(cos_theta, l, 0);
    else
      return std::sqrt(2) * N_lm(l, m) * evaluate_associated_legendre_polynomial(cos_theta, l, m) * std::cos(m * phi);
  }
}; // class RealSphericalHarmonicsBase<DomainFieldType, 3, ...>

template <class DomainFieldType, class RangeFieldType, size_t order>
class FullSphericalHarmonics<DomainFieldType, 3, RangeFieldType, order>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, (order + 1) * (order + 1), 1>
{
  static const size_t dimDomain = 3;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  using BaseType::DomainType;
  using BaseType::RangeType;

  virtual RangeType evaluate(const DomainType& v) const override
  {
    const auto v_spherical = CoordinateConverter::to_spherical(v);
    return evaluate_in_spherical_coords(v_spherical);
  } // ... evaluate(...)

  RangeType evaluate_in_spherical_coords(const FieldVector<DomainFieldType, 2>& coords) const
  {
    const DomainFieldType theta = coords[0];
    const DomainFieldType phi = coords[1];
    RangeType ret(0);
    for (int ll = 0; ll <= order; ++ll)
      for (int mm = -ll; mm <= ll; ++mm)
        ret[pos(ll, mm)] = boost::math::spherical_harmonic(ll, mm, theta, phi);
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

  FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimDomain> ret(MatrixType(0));
    const auto& Bx = ret[0];
    const auto& By = ret[1];
    const auto& Bz = ret[2];
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -l1; std::abs(m1) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -l2; std::abs(m2) <= l2; ++m2) {
            if (l1 == l2 + 1 && m1 == m2 + 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::complex(0, 1) * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            }
            if (l1 == l2 - 1 && m1 == m2 + 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::complex(0, 1) * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            }
            if (l1 == l2 + 1 && m1 == m2 - 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::complex(0, 1) * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            }
            if (l1 == l2 - 1 && m1 == m2 - 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::complex(0, 1) * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            }
            if (m1 == m2 && l1 == l2 + 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2 + 1) * (l2 + m2 + 1) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (m1 == m2 && l1 == l2 - 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2) * (l2 + m2) / ((2 * l2 - 1) * (2 * l2 + 1)));
          } // m2
        } // l2
      } // m1
    } // l1
    return ret;
  }

private:
  // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
  // Each l has 2l+1 values of m, so (l, m) has position
  // (\sum_{k=0}^{l-1} (2k+1)) + (m+l) = l^2 + m + l
  size_t pos(const size_t l, const int m) const
  {
    return l * l + m + l;
  }
}; // class FullRealSphericalHarmonics<DomainFieldType, 3, ...>

template <class DomainFieldType, class RangeFieldType, size_t order>
class PositiveSphericalHarmonics<DomainFieldType, 3, RangeFieldType, order>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, (order + 1) * (order + 2) / 2, 1>
{
  static const size_t dimDomain = 3;
  static const size_t dimRange = (order + 1) * (order + 2) / 2;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  using BaseType::DomainType;
  using BaseType::RangeType;

  virtual RangeType evaluate(const DomainType& v) const override
  {
    const auto v_spherical = CoordinateConverter::to_spherical(v);
    return evaluate_in_spherical_coords(v_spherical);
  } // ... evaluate(...)

  RangeType evaluate_in_spherical_coords(const FieldVector<DomainFieldType, 2>& coords) const
  {
    const DomainFieldType theta = coords[0];
    const DomainFieldType phi = coords[1];
    RangeType ret(0);
    for (size_t ll = 0; ll <= order; ++ll)
      for (int mm = 0; size_t(mm) <= ll; ++mm)
        ret[pos(ll, mm)] = boost::math::spherical_harmonic(ll, mm, theta, phi);
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

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimDomain> ret(MatrixType(0));
    const auto& Bx = ret[0];
    const auto& By = ret[1];
    const auto& Bz = ret[2];
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (size_t m1 = 0; m1 <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (size_t m2 = 0; m2 <= l2; ++m2) {
            if (l1 == l2 + 1 && m1 == m2 + 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::complex(0, 1) * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            }
            if (l1 == l2 - 1 && m1 == m2 + 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::complex(0, 1) * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            }
            if (l1 == l2 + 1 && m1 == m2 - 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::complex(0, 1) * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            }
            if (l1 == l2 - 1 && m1 == m2 - 1) {
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              By[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::complex(0, 1) * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            }
            if (m1 == m2 && l1 == l2 + 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2 + 1) * (l2 + m2 + 1) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (m1 == m2 && l1 == l2 - 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2) * (l2 + m2) / ((2 * l2 - 1) * (2 * l2 + 1)));
          } // m2
        } // l2
      } // m1
    } // l1
    return ret;
  }

private:
  // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
  // Each l has l+1 non-negative values of m, so (l, m) has position
  // (\sum_{k=0}^{l-1} (l+1)) + m = l(l+1)/2 + m
  size_t pos(const size_t l, const int m) const
  {
    return l * (l + 1) / 2 + m;
  }
}; // class PositiveSphericalHarmonics<DomainFieldType, 3, ...>


template <class DomainFieldType, class RangeFieldType, size_t order>
class FullRealSphericalHarmonics<DomainFieldType, 3, RangeFieldType, order>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, (order + 1) * (order + 1), 1>,
      public RealSphericalHarmonicsBase<DomainFieldType, 3, RangeFieldType>
{
  static const size_t dimDomain = 3;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;
  typedef RealSphericalHarmonicsBase<DomainFieldType, 3, RangeFieldType> SecondBaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using BaseType::DomainType;
  using BaseType::RangeType;

  virtual RangeType evaluate(const DomainType& v) const override
  {
    const auto v_spherical = CoordinateConverter::to_spherical(v);
    return evaluate_in_spherical_coords(v_spherical);
  } // ... evaluate(...)

  RangeType evaluate_in_spherical_coords(const FieldVector<DomainFieldType, 2>& coords) const
  {
    const DomainFieldType theta = coords[0];
    const DomainFieldType phi = coords[1];
    RangeType ret(0);
    for (int ll = 0; ll <= order; ++ll)
      for (int mm = -ll; mm <= ll; ++mm)
        ret[pos(ll, mm)] = evaluate_real_spherical_harmonics(theta, phi, ll, mm);
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

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimDomain> ret(MatrixType(0));
    const auto& Bx = ret[0];
    const auto& By = ret[1];
    const auto& Bz = ret[2];
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -l1; std::abs(m1) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -l2; std::abs(m2) <= l2; ++m2) {
            if (l1 == l2 + 1 && m1 == m2 + 1)
              Bx[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (l1 == l2 - 1 && m1 == m2 + 1)
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            if (l1 == l2 + 1 && m1 == m2 - 1)
              Bx[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (l1 == l2 - 1 && m1 == m2 - 1)
              Bx[pos(l1, m1)][pos(l2, m2)] = 0.5 * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            if (l1 == l2 + 1 && m1 == -m2 - 1)
              By[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * sgn(m2) * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (l1 == l2 - 1 && m1 == -m2 - 1)
              By[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * sgn(m2) * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            if (l1 == l2 + 1 && m1 == -m2 + 1)
              By[pos(l1, m1)][pos(l2, m2)] =
                  0.5 * sgn(m2) * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (l1 == l2 - 1 && m1 == -m2 + 1)
              By[pos(l1, m1)][pos(l2, m2)] =
                  -0.5 * sgn(m2) * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
            if (m1 == m2 && l1 == l2 + 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2 + 1) * (l2 + m2 + 1) / ((2 * l2 + 1) * (2 * l2 + 3)));
            if (m1 == m2 && l1 == l2 - 1)
              Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2) * (l2 + m2) / ((2 * l2 - 1) * (2 * l2 + 1)));
          } // m2
        } // l2
      } // m1
    } // l1
    return ret;
  }

private:
  template <typename T>
  int sgn(T val)
  {
    return (T(0) < val) - (val < T(0));
  }
  // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
  // Each l has 2l+1 values of m, so (l, m) has position
  // (\sum_{k=0}^{l-1} (2k+1)) + (m+l) = l^2 + m + l
  size_t pos(const int l, const int m) const
  {
    return size_t(l * l + m + l);
  }

  using SecondBaseType::evaluate_real_spherical_harmonics;
}; // class FullRealSphericalHarmonics<DomainFieldType, 3, ...>

template <class DomainFieldType, class RangeFieldType, size_t order>
class EvenRealSphericalHarmonics<DomainFieldType, 3, RangeFieldType, order>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, (order + 1) * (order + 2) / 2, 1>,
      public RealSphericalHarmonicsBase<DomainFieldType, 3, RangeFieldType>
{
  static const size_t dimDomain = 3;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;
  typedef RealSphericalHarmonicsBase<DomainFieldType, 3, RangeFieldType> SecondBaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using BaseType::DomainType;
  using BaseType::RangeType;

  virtual RangeType evaluate(const DomainType& v) const override
  {
    const auto v_spherical = CoordinateConverter::to_spherical(v);
    return evaluate_in_spherical_coords(v_spherical);
  } // ... evaluate(...)

  RangeType evaluate_in_spherical_coords(const FieldVector<DomainFieldType, 2>& coords) const
  {
    const DomainFieldType theta = coords[0];
    const DomainFieldType phi = coords[1];
    RangeType ret(0);
    for (int ll = 0; ll <= order; ++ll)
      for (int mm = -ll; mm <= ll; ++mm)
        if (!((mm + ll) % 2))
          ret[pos(ll, mm)] = evaluate_real_spherical_harmonics(theta, phi, ll, mm);
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

  MatrixType mass_matrix_inverse() const
  {
    return mass_matrix();
  }

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimDomain> ret(MatrixType(0));
    const auto& Bx = ret[0];
    const auto& By = ret[1];
    const auto& Bz = ret[2];
    for (size_t l1 = 0; l1 <= order; ++l1) {
      for (int m1 = -l1; std::abs(m1) <= l1; ++m1) {
        for (size_t l2 = 0; l2 <= order; ++l2) {
          for (int m2 = -l2; std::abs(m2) <= l2; ++m2) {
            if (!((m + l) % 2) && !((m2 + l2) % 2)) {
              if (l1 == l2 + 1 && m1 == m2 + 1)
                Bx[pos(l1, m1)][pos(l2, m2)] =
                    -0.5 * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              if (l1 == l2 - 1 && m1 == m2 + 1)
                Bx[pos(l1, m1)][pos(l2, m2)] =
                    0.5 * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              if (l1 == l2 + 1 && m1 == m2 - 1)
                Bx[pos(l1, m1)][pos(l2, m2)] =
                    0.5 * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              if (l1 == l2 - 1 && m1 == m2 - 1)
                Bx[pos(l1, m1)][pos(l2, m2)] =
                    0.5 * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              if (l1 == l2 + 1 && m1 == -m2 - 1)
                By[pos(l1, m1)][pos(l2, m2)] =
                    0.5 * sgn(m2) * std::sqrt((l2 + m2 + 1) * (l2 + m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              if (l1 == l2 - 1 && m1 == -m2 - 1)
                By[pos(l1, m1)][pos(l2, m2)] =
                    -0.5 * sgn(m2) * std::sqrt((l2 - m2) * (l2 - m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              if (l1 == l2 + 1 && m1 == -m2 + 1)
                By[pos(l1, m1)][pos(l2, m2)] =
                    0.5 * sgn(m2) * std::sqrt((l2 - m2 + 1) * (l2 - m2 + 2) / ((2 * l2 + 1) * (2 * l2 + 3)));
              if (l1 == l2 - 1 && m1 == -m2 + 1)
                By[pos(l1, m1)][pos(l2, m2)] =
                    -0.5 * sgn(m2) * std::sqrt((l2 + m2) * (l2 + m2 - 1) / ((2 * l2 - 1) * (2 * l2 + 1)));
              if (m1 == m2 && l1 == l2 + 1)
                Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2 + 1) * (l2 + m2 + 1) / ((2 * l2 + 1) * (2 * l2 + 3)));
              if (m1 == m2 && l1 == l2 - 1)
                Bz[pos(l1, m1)][pos(l2, m2)] = std::sqrt((l2 - m2) * (l2 + m2) / ((2 * l2 - 1) * (2 * l2 + 1)));
            }
          } // m2
        } // l2
      } // m1
    } // l1
    return ret;
  }

private:
  // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
  // Each l has l+1 values of m (as only m s.t. m+l is even are considered), so (l, m) has position
  // (\sum_{k=0}^{l-1} (k+1)) + (m+l)/2 = l(l+1)/2 + (l+m)/2
  size_t pos(const int l, const int m) const
  {
    return size_t(l * (l + 1) / 2 + (m + l) / 2);
  }

  using SecondBaseType::evaluate_real_spherical_harmonics;
}; // class EvenRealSphericalHarmonics<DomainFieldType, 3, ...>


template <class DomainFieldType, class RangeFieldType, size_t dimRange, size_t dimRangeCols>
class HatFunctions<DomainFieldType, 1, RangeFieldType, dimRange, dimRangeCols>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, dimRange, dimRangeCols>
{
  static const size_t dimDomain = 1;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using BaseType::DomainType;
  using BaseType::RangeType;
  typedef RangeType TriangulationType;

  HatFunctions(const TriangulationType& triangulation, const QuadratureType& /*quadrature*/ = QuadratureType())
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
      if (ii < dimRange - 1 && XT::Common::FloatCmp::ge(v, triangulation_[ii])
          && XT::Common::FloatCmp::le(v, triangulation_[ii + 1]))
        ret[ii] = (v - triangulation_[ii + 1]) / (triangulation_[ii] - triangulation_[ii + 1]);
      if (ii > 0 && XT::Common::FloatCmp::ge(v, triangulation_[ii - 1])
          && XT::Common::FloatCmp::le(v, triangulation_[ii]))
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
    return tridiagonal_matrix_inverse(mass_matrix);
  }

  // returns matrix with entries <v h_i h_j>
  virtual MatrixType mass_matrix_with_v() const override
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

private:
  const TriangulationType& triangulation_;
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

template <class DomainFieldType, class RangeFieldType, size_t dimRange, size_t dimRangeCols>
class HatFunctions<DomainFieldType, 3, RangeFieldType, dimRange, dimRangeCols>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, dimRange, dimRangeCols>
{
  static const size_t dimDomain = 3;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef typename SphericalTriangulation<DomainFieldType> TriangulationType;
  using BaseType::DomainType;
  using BaseType::RangeType;

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
    for (const auto& face : triangulation.faces()) {
      const auto& vertices = face->vertices();
      DomainType barycentric_coords(0);
      success = calculate_barycentric_coordinates(v, vertices, barycentric_coords);
      if (success) {
        for (size_t ii = 0; ii < 3; ++ii) {
          ret[vertices[ii]->index()] = barycentric_coords[ii];
        }
        break;
      }
    } // faces
    assert(success);
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override
  {
    RangeType ret(0);
    for (const auto& quad_point : quadrature_) {
      const auto v = quad_point.position();
      const auto basis_evaluated = evaluate(v);
      const auto weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        ret[nn] += basis_evaluated[nn] * weight;
    } // quadrature
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType A;
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
    return tridiagonal_matrix_inverse(mass_matrix());
  }

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimDomain> B(MatrixType(0));
    for (const auto& quad_point : quadrature_) {
      const auto& v = quad_point.position();
      const auto basis_evaluated = evaluate(v);
      const auto& weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        for (size_t mm = 0; mm < dimRange; ++mm)
          for (size_t dd = 0; dd < dimDomain; ++dd)
            B[dd][nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * v[dd] * weight;
    } // quadrature
    return B;
  } // ... create_flux_config(...)

protected:
  template <class VertexVectorType>
  bool calculate_barycentric_coordinates(const DomainType& v, const VertexVectorType& vertices, DomainType& ret)
  {
    Dune::FieldMatrix<double, 3, 3> gradients(0);
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
      gradients[ii] *= std::acos(scalar_prod) / std::sqrt(1. - std::pow(scalar_prod, 2));
    } // ii
    // calculate barycentric coordinates for 0 w.r.t the points g_i
    const auto& g0 = gradients[0];
    const auto& g1 = gradients[1];
    const auto& g2 = gradients[2];
    auto g0_minus_g2 = g0;
    auto g1_minus_g2 = g1;
    g0_minus_g2 -= g2;
    g1_minus_g2 -= g2;
    Dune::FieldMatrix<double, 2, 2> A;
    Dune::FieldVector<double, 2> solution;
    Dune::FieldVector<double, 2> rhs;
    // (ii, jj) = (0, 1), (0, 2), (1, 2)
    for (size_t ii = 0; ii < 2; ++ii) {
      for (size_t jj = ii + 1; jj < 3; ++jj) {
        A[0][0] = g0_minus_g2[ii];
        A[1][0] = g0_minus_g2[jj];
        A[0][1] = g1_minus_g2[ii];
        A[1][1] = g1_minus_g2[jj];
        double det = A.determinant();
        if (XT::Common::FloatCmp::eq(det, 0.))
          break;
        rhs[0] = -g2[ii];
        rhs[1] = -g2[jj];
        A.solve(solution, rhs);
        if (XT::Common::FloatCmp::lt(solution[0], 0.) || XT::Common::FloatCmp::lt(solution[1], 0.))
          return false;
        ret[0] = solution[0];
        ret[1] = solution[1];
        ret[2] = 1. - ret[0] - ret[1];
        if (XT::Common::FloatCmp::lt(ret[2], 0.))
          return false;
        return true;
      }
    }
    return false;
  } // bool calculate_barycentric_coordinates(...)

  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
}; // class HatFunctions<DomainFieldType, 3, ...>


template <class DomainFieldType, class RangeFieldType, size_t dimRange, size_t dimRangeCols>
class PiecewiseMonomials<DomainFieldType, 1, RangeFieldType, dimRange, dimRangeCols, 1>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, dimRange, dimRangeCols>
{
  static_assert(!(dimRange % 2), "dimRange has to be even!");
  static const size_t dimDomain = 1;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef FieldVector<DomainFieldType, dimRange + 1> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  PiecewiseMonomials(const TriangulationType& triangulation, const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    RangeType ret;
    for (size_t ii = 0; ii < dimRange / 2 + 1; ++ii)
      ret[ii] = -1. + 4. * ii / dimRange;
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange; ii += 2) {
      if (XT::Common::FloatCmp::ge(v, triangulation_[ii / 2])
          && XT::Common::FloatCmp::le(v, triangulation_[ii / 2 + 1])) {
        ret[ii] = 1;
        ret[ii + 1] = v;
      }
    }
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange; ii += 2) {
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
      M[2 * ii + 1][2 * ii] = ret[2 * ii][2 * ii + 1];
    }
    return ret;
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    MatrixType B(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      B[2 * ii][2 * ii] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      B[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 4) - std::pow(triangulation_[ii], 4)) / 4.;
      B[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      B[2 * ii + 1][2 * ii] = ret[2 * ii][2 * ii + 1];
    }
    return FieldVector<MatrixType, dimDomain>(B);
  }

private:
  using BaseType::triangulation_;
}; // class PiecewiseMonomials<DomainFieldType, 1, ...>


template <class DomainFieldType, class RangeFieldType, size_t dimRange, size_t dimRangeCols>
class PiecewiseMonomials<DomainFieldType, 3, RangeFieldType, dimRange, dimRangeCols>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, dimRange, dimRangeCols>
{
  static const size_t dimDomain = 3;
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef typename SphericalTriangulation<DomainFieldType> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

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
    FieldMatrix<double, 3, 3> vertices_matrix;
    FieldMatrix<double, 3, 3> determinant_matrix;
    for (const auto& face : triangulation_.faces()) {
      // vertices are ordered counterclockwise, so if the points is inside the spherical triangle,
      // the coordinate system formed by two adjacent vertices and v is always right-handed, i.e.
      // the triple product is positive
      const auto& vertices = face->vertices();
      for (size_t ii = 0; ii < 3; ++ii)
        vertices_matrix[ii] = vertices[ii]->position();
      bool v_in_this_facet = true;
      // the triple products that need to be positive are the determinants of the matrices (v1, v2, v), (v2, v3, v),
      // (v3,
      // v1, v), where vi is the ith
      // vertex. Swapping two columns changes the sign of det, the matrices used below all have an even number of column
      // swaps
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
    RangeType ret(0);
    for (const auto& quad_point : quadrature_) {
      const auto basis_evaluated = evaluate(quad_point.position());
      for (size_t nn = 0; nn < dimRange; ++nn)
        ret[nn] += basis_evaluated[nn] * quad_point.weight();
    } // quadrature
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(0);
    for (const auto& quad_point : quadrature_) {
      const auto basis_evaluated = evaluate(quad_point.weight());
      for (size_t nn = 0; nn < dimRange; ++nn)
        for (size_t mm = 0; mm < dimRange; ++mm)
          M[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * weight;
    } // quadrature
    return M;
  } // ... mass_matrix()

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimDomain> B(MatrixType(0));
    for (const auto& quad_point : quadrature_) {
      const auto v = quad_point.position();
      const auto basis_evaluated = evaluate(v);
      const auto weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        for (size_t mm = 0; mm < dimRange; ++mm)
          for (size_t dd = 0; dd < dimDomain; ++dd)
            B[dd][nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * v[dd] * weight;
    } // quadrature
    return B;
  }

private:
  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
}; // class PiecewiseMonomials<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH
