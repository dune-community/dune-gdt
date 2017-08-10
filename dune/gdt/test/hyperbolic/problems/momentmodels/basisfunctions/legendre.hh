// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

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
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

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
    MatrixType M(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 2. / (2. * rr + 1.);
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    MatrixType Minv(dimRange, dimRange, 0);
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
          B[rr][cc] = 2 * rr / (4. * rr * rr - 1.);
        else if (cc == rr + 1)
          B[rr][cc] = (2 * rr + 2.) / ((2. * rr + 1.) * (2. * rr + 3));
      }
    }
    return B;
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


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEGENDRE_HH
