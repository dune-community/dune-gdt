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


// TODO: use complex arithmetic, currently only usable for Pn Models in 2D, test for only_positive = false
template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim = 3, bool only_positive = true>
class SphericalHarmonicsMomentBasis
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
  using typename BaseType::QuadraturesType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;
  static_assert(order <= std::numeric_limits<int>::max(), "");

  SphericalHarmonicsMomentBasis(const QuadraturesType& quadratures)
    : BaseType(quadratures)
  {
  }

  SphericalHarmonicsMomentBasis(const size_t quad_order = 2 * order + 2,
                                const size_t DXTC_DEBUG_ONLY(quad_refinements) = 0)
    : BaseType(OctantQuadrature<DomainFieldType>::get(quad_order))
  {
    assert(quad_refinements == 0 && "Refinement of the quadrature intervals not implemented for this basis!");
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
    // TODO: use complex arithmetic, remove real() call
    for (unsigned int ll = 0; ll <= static_cast<unsigned int>(order); ++ll)
      for (int mm = only_positive ? 0 : -static_cast<int>(ll); mm <= static_cast<int>(ll); ++mm)
        ret[helper<only_positive>::pos(ll, mm)] = boost::math::spherical_harmonic(ll, mm, theta, phi).real();
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated(const bool /*use_fine_quadratures*/ = false) const override
  {
    RangeType ret(0);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  virtual MatrixType mass_matrix(const bool /*use_fine_quadratures*/ = false) const override
  {
    MatrixType M(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 1;
    return M;
  }

  virtual MatrixType mass_matrix_inverse(const bool /*use_fine_quadratures*/ = false) const override
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

  virtual RangeType alpha_iso() const override final
  {
    RangeType ret(0.);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    return u[0] * std::sqrt(4 * M_PI);
  }

  virtual std::string short_id() const override final
  {
    return "shm";
  }

private:
  static RangeFieldType A_lm(const int l, const int m)
  {
    assert(std::abs(m) <= l);
    return std::sqrt((l + m) * (l - m) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static RangeFieldType B_lm(const int l, const int m)
  {
    assert(std::abs(m) <= l);
    return std::sqrt((l + m) * (l + m - 1.) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static MatrixType create_Bx()
  {
    MatrixType Bx(dimRange, dimRange, 0);
    const auto& pos = helper<only_positive>::pos;
    for (int l1 = 0; l1 <= static_cast<int>(order); ++l1) {
      for (int m1 = only_positive ? 0 : -int(l1); std::abs(m1) <= l1; ++m1) {
        for (int l2 = 0; l2 <= static_cast<int>(order); ++l2) {
          for (int m2 = -int(l2); std::abs(m2) <= l2; ++m2) {
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
    for (int l1 = 0; l1 <= static_cast<int>(order); ++l1) {
      for (int m1 = only_positive ? 0. : -l1; std::abs(m1) <= l1; ++m1) {
        for (int l2 = 0; l2 <= static_cast<int>(order); ++l2) {
          size_t row = pos(l1, m1);
          if (l1 == l2 + 1 && std::abs(m1) < l1) {
            size_t col = pos(l2, m1); // m1 == m2, else matrix entry is 0
            Bz[row][col] += A_lm(l2 + 1, m1);
          }
          if (l1 == l2 - 1) {
            size_t col = pos(l2, m1); // m1 == m2, else matrix entry is 0
            Bz[row][col] += A_lm(l2, m1);
          }
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
      const int ret = (l * l + m + l);
      assert(ret >= 0 && std::abs(m) <= l);
      return static_cast<size_t>(ret);
    }
  };

  template <class anything>
  struct helper<true, anything>
  {
    // Converts a pair (l, m) to a vector index. The vector is ordered by l first, then by m.
    // Each l has l+1 non-negative values of m, so (l, m) has position
    // (\sum_{k=0}^{l-1} (l+1)) + m = l(l+1)/2 + m
    static size_t pos(const int l, const int m)
    {
      const int ret = l * (l + 1) / 2 + m;
      assert(std::abs(m) <= l && ret >= 0);
      return static_cast<size_t>(ret);
    }
  };
}; // class SphericalHarmonicsMomentBasis<DomainFieldType, 3, ...>


template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim = 3, bool only_even = false>
class RealSphericalHarmonicsMomentBasis
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
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  RealSphericalHarmonicsMomentBasis(const QuadraturesType& quadratures)
    : BaseType(quadratures)
  {
  }

  RealSphericalHarmonicsMomentBasis(const size_t quad_order = 2 * order + 2,
                                    const size_t DXTC_DEBUG_ONLY(quad_refinements) = 0)
    : BaseType(OctantQuadrature<DomainFieldType>::get(quad_order))
  {
    assert(quad_refinements == 0 && "Refinement of the quadrature intervals not implemented for this basis!");
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
    for (int ll = 0; ll <= static_cast<int>(order); ++ll)
      for (int mm = -ll; mm <= ll; ++mm)
        if (!only_even || !((mm + ll) % 2))
          ret[helper<only_even>::pos(ll, mm)] = evaluate_lm(theta, phi, ll, mm);
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated(const bool /*use_fine_quadratures*/ = false) const override
  {
    RangeType ret(0.);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  virtual MatrixType mass_matrix(const bool /*use_fine_quadratures*/ = false) const override
  {
    MatrixType M(dimRange, dimRange, 0.);
    for (size_t rr = 0; rr < dimRange; ++rr)
      M[rr][rr] = 1;
    return M;
  }

  virtual MatrixType mass_matrix_inverse(const bool /*use_fine_quadratures*/ = false) const override
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

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      component_visualizer(0, u_n, filename_prefix, ii, std::sqrt(4 * M_PI));
    };
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) { return XT::Common::to_string(val[0] * std::sqrt(4 * M_PI), 15); };
  } // ... stringifier()

  virtual RangeType alpha_iso() const override final
  {
    RangeType ret(0.);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    return u[0] * std::sqrt(4 * M_PI);
  }

  virtual std::string short_id() const override final
  {
    return "rhm";
  }

  RangeType integrate_dirac_at(const DomainType& dirac_position) const
  {
    return evaluate(dirac_position);
  }

private:
  static RangeFieldType A_lm(const int l, const int m)
  {
    assert(std::abs(m) <= l);
    return std::sqrt((l + m) * (l - m) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static RangeFieldType B_lm(const int l, const int m)
  {
    assert(std::abs(m) <= l);
    return std::sqrt((l + m) * (l + m - 1.) / ((2. * l + 1.) * (2. * l - 1.)));
  }

  static MatrixType create_Bx()
  {
    MatrixType Bx(dimRange, dimRange, 0.);
    const auto& pos = helper<only_even>::pos;
    for (int l1 = 0; l1 <= static_cast<int>(order); ++l1) {
      for (int m1 = -int(l1); std::abs(m1) <= l1; ++m1) {
        for (int l2 = 0; l2 <= static_cast<int>(order); ++l2) {
          for (int m2 = -int(l2); std::abs(m2) <= l2; ++m2) {
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
    for (int l1 = 0; l1 <= static_cast<int>(order); ++l1) {
      for (int m1 = -int(l1); std::abs(m1) <= l1; ++m1) {
        for (int l2 = 0; l2 <= static_cast<int>(order); ++l2) {
          for (int m2 = -int(l2); std::abs(m2) <= l2; ++m2) {
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
    for (int l1 = 0; l1 <= static_cast<int>(order); ++l1) {
      for (int m1 = -l1; std::abs(m1) <= l1; ++m1) {
        for (int l2 = 0; l2 <= static_cast<int>(order); ++l2) {
          for (int m2 = -l2; std::abs(m2) <= l2; ++m2) {
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
    static size_t pos(const int l, const int m)
    {
      const int ret = l * l + m + l;
      assert(ret >= 0 && std::abs(m) <= l);
      return static_cast<size_t>(ret);
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
      const int ret = l * (l + 1) / 2 + (m + l) / 2;
      assert(std::abs(m) <= l && ret >= 0);
      return static_cast<size_t>(ret);
    }
  };

  // Notation from Garrett, Hauck, "A Comparison of Moment Closures for Linear Kinetic Transport Equations: The Line
  // Source Benchmark",
  // http://www.tandfonline.com/doi/full/10.1080/00411450.2014.910226?src=recsys&, Section 4.1
  RangeFieldType N_lm(const int l, const int m) const
  {
    static constexpr auto frac_4pi = 1. / (4. * M_PI);
    assert(l >= 0 && m >= 0 && m <= l);
    // return std::sqrt((2. * l + 1.) * XT::Common::factorial(l - m) / (XT::Common::factorial(l + m) * 4. * M_PI));
    auto factor = 1.;
    for (int ii = l - m + 1; ii <= l + m; ++ii)
      factor *= ii;
    factor = 1. / factor;
    return std::sqrt((2. * l + 1.) * factor * frac_4pi);
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

  using BaseType::quadratures_;
}; // class RealSphericalHarmonicsMomentBasis<DomainFieldType, 3, ...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_SPHERICALHARMONICS_HH
