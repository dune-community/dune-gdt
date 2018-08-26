// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEGENDRE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEGENDRE_HH

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class DomainFieldType, class RangeFieldType, size_t order, size_t dimRangeCols = 1>
class LegendrePolynomials : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, order + 1, dimRangeCols>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = order + 1;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;
  using TriangulationType = typename BaseType::Triangulation1dType;

  LegendrePolynomials(const QuadraturesType& quadratures)
    : BaseType(quadratures)
  {
  }

  LegendrePolynomials(const size_t quad_order = 31, const size_t quad_refinements = 0)
    : BaseType(BaseType::gauss_lobatto_quadratures(std::pow(2, quad_refinements), quad_order))
  {
  }

  static std::string static_id()
  {
    return "legendre";
  }

  virtual RangeType evaluate(const DomainType& v) const override
  {
    RangeType ret;
    ret[0] = 1.;
    if (dimRange > 1)
      ret[1] = v[0];
    for (size_t ii = 2; ii < dimRange; ++ii)
      ret[ii] = ((2. * ii - 1.) * v[0] * ret[ii - 1] - (ii - 1.) * ret[ii - 2]) / ii;
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
    MatrixType M(dimRange, dimRange, 0.);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 2. / (2. * rr + 1.);
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    MatrixType Minv(dimRange, dimRange, 0.);
    for (size_t rr = 0; rr < dimRange; ++rr)
      Minv[rr][rr] = (2. * rr + 1.) / 2.;
    return Minv;
  }

  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    MatrixType B(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr) {
      for (size_t cc = 0; cc < dimRange; ++cc) {
        if (cc == rr - 1)
          B[rr][cc] = 2. * rr / (4. * rr * rr - 1.);
        else if (cc == rr + 1)
          B[rr][cc] = (2. * rr + 2.) / ((2. * rr + 1.) * (2. * rr + 3));
      }
    }
    return B;
  }

  // returns matrices with entries <v h_i h_j>_- and <v h_i h_j>_+
  virtual FieldVector<FieldVector<MatrixType, 2>, 1> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, 1> ret(FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    auto mm_with_v = mass_matrix_with_v();
    auto& ret_neg = ret[0][0];
    auto& ret_pos = ret[0][1];
    MatrixType mass_matrix_pos(dimRange + 1, dimRange + 1, 0.); // we need to calculate up to P_N
    size_t N = dimRange;
    // calculate <P_n P_m>_+ first
    for (size_t nn = 0; nn <= N; ++nn) {
      for (size_t mm = 0; mm <= N; ++mm) {
        if ((nn + mm) % 2) { // nn and mm have different parity
          if (nn % 2) // n odd
            mass_matrix_pos[nn][mm] = fmn(static_cast<int>(mm), static_cast<int>(nn));
          else // n even
            mass_matrix_pos[nn][mm] = fmn(static_cast<int>(nn), static_cast<int>(mm));
        } else { // nn and mm are both odd or both even
          mass_matrix_pos[nn][mm] = (nn == mm) ? 1. / (2. * nn + 1.) : 0.;
        }
      } // mm
    } // nn
    // now calculate <v P_n P_m>_+ and <v P_n P_m>_-
    for (size_t nn = 0; nn < N; ++nn) {
      for (size_t mm = 0; mm < N; ++mm) {
        ret_pos[nn][mm] = (nn + 1.) / (2. * nn + 1.) * mass_matrix_pos[nn + 1][mm]
                          + ((nn > 0) ? nn / (2. * nn + 1.) * mass_matrix_pos[nn - 1][mm] : 0.);
        ret_neg[nn][mm] = mm_with_v[0][nn][mm] - ret_pos[nn][mm];
      } // mm
    } // nn
    return ret;
  }

  virtual MatrixType reflection_matrix(const DomainType& n) const override final
  {
    MatrixType ret(dimRange, dimRange, 0);
    for (size_t ii = 0; ii < dimDomain; ++ii)
      if (XT::Common::FloatCmp::ne(n[ii], 0.))
        if (XT::Common::FloatCmp::ne(std::abs(n[ii]), 1.))
          DUNE_THROW(NotImplemented, "Implemented only for +-e_i where e_i is the i-th canonical basis vector!");
    const auto mass_mat = mass_matrix();
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        ret[ii][jj] = std::pow(-1, ii) * mass_mat[ii][jj];
    ret.rightmultiply(mass_matrix_inverse());
    return ret;
  }

  MatrixType S() const
  {
    MatrixType S(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      S[rr][rr] = -2. * rr * (rr + 1.) / (2 * rr + 1);
    return S;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      component_visualizer<DiscreteFunctionType, dimRange, 0>(u_n, filename_prefix, ii);
    };
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) { return XT::Common::to_string(val[0], 15); };
  } // ... stringifier()

  virtual RangeType alpha_iso() const override final
  {
    RangeType ret(0.);
    ret[0] = 1.;
    return ret;
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    return u[0];
  }

  virtual std::string short_id() const override final
  {
    return "leg";
  }

private:
  static RangeFieldType fmn(const int m, const int n)
  {
    assert(!(m % 2));
    assert(n % 2);
    // The factorials overflow for large m or n, so do it more complicated
    // return (std::pow(-1., (m + n + 1) / 2) * XT::Common::factorial(m) * XT::Common::factorial(n))
    //       / (std::pow(2., m + n - 1) * (m - n) * (m + n + 1) * std::pow(XT::Common::factorial(m / 2), 2)
    //          * std::pow(XT::Common::factorial((n - 1) / 2), 2));

    RangeFieldType ret = 1;
    int m_factor = m;
    int n_factor = n;
    int m_divisor = m / 2;
    int n_divisor = (n - 1) / 2;
    size_t max_mn = static_cast<size_t>(std::max(m, n));
    FieldVector<std::vector<RangeFieldType>, 4> factors((std::vector<RangeFieldType>(max_mn)));
    assert(std::max(m, n) >= 0);
    for (size_t ii = 0; ii < max_mn; ++ii) {
      factors[0][ii] = m_factor > 0 ? m_factor-- : 1.;
      factors[1][ii] = n_factor > 0 ? n_factor-- : 1.;
      factors[2][ii] = m_divisor > 0 ? 1. / std::pow(m_divisor--, 2) : 1.;
      factors[3][ii] = n_divisor > 0 ? 1. / std::pow(n_divisor--, 2) : 1.;
    }
    for (size_t ii = 0; ii < max_mn; ++ii) {
      ret *= factors[0][ii] * factors[1][ii] * factors[2][ii] * factors[3][ii] / 2.;
    }
    ret *= std::pow(-1., (m + n + 1) / 2) / ((m - n) * (m + n + 1) * std::pow(2., m + n - 1 - std::max(m, n)));
    return ret;
  } // ... fmn(...)
}; // class LegendrePolynomials<DomainFieldType, 1, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEGENDRE_HH
