// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/math.hh>

#include <dune/xt/la/eigen-solver.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class SourceBeamPn : public KineticTransportEquationBase<E, MomentBasisImp>
{
  using BaseType = KineticTransportEquationBase<E, MomentBasisImp>;
  using ThisType = SourceBeamPn;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::BoundaryDistributionType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::FluxType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::GenericScalarFunctionType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::ScalarFunctionType;

  SourceBeamPn(const MomentBasis& basis_functions,
               const RangeFieldType psi_vac = 5e-7,
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const bool is_mn_model = false)
    : BaseType(basis_functions, psi_vac, grid_cfg)
    , is_mn_model_(is_mn_model)
  {}

  static std::string static_id()
  {
    return "sourcebeampn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[3.0]";
    grid_config["num_elements"] = "[300]";
    grid_config["overlap_size"] = "[2]";
    return grid_config;
  }

  // Boundary value of kinetic equation is \frac{g}{<g>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, with g(v) = exp(-10^5(v-1)^2), so n-th component of boundary value has to be
  // \frac{<base_n(v)*g(v)>}{<g>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  std::unique_ptr<BoundaryValueType> boundary_values() const override final
  {
    return std::make_unique<GenericFunctionType>(
        1, [&](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) {
          if (x[0] < 1.5) {
            helper<MomentBasis>::get_left_boundary_values(basis_functions_, psi_vac_, is_mn_model_, ret);
          } else {
            ret = basis_functions_.integrated();
            ret *= psi_vac_;
          }
        });
  } // ... boundary_values()

  BoundaryDistributionType boundary_distribution() const override final
  {
    return [this](const DomainType& x) -> std::function<RangeFieldType(const DomainType&)> {
      if (x[0] > 1.5)
        return [this](const DomainType& /*v*/) { return this->psi_vac_; };
      else
        return [this](const DomainType& v) {
          const RangeFieldType val = std::exp(-1e5 * std::pow(v[0] - 1., 2)) / helper_base::denominator();
          return val + psi_vac_;
        };
    };
  }

  RangeReturnType
  kinetic_boundary_flux(const DomainType& x, const RangeFieldType& n, const size_t dd) const override final
  {
    return helper<MomentBasis>::get_kinetic_boundary_flux(basis_functions_, psi_vac_, is_mn_model_, x, n, dd, *this);
  }

  RangeReturnType left_boundary_value() const
  {
    DynamicRangeType ret;
    helper<MomentBasis>::get_left_boundary_values(basis_functions_, psi_vac_, is_mn_model_, ret);
    return ret;
  }

  RangeFieldType t_end() const override
  {
    return 2.5;
  }

  // sigma_a = 1 if x <= 2, 0 else
  std::unique_ptr<ScalarFunctionType> sigma_a() const override
  {
    return std::make_unique<GenericScalarFunctionType>(
        0, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] > 2 ? 0. : 1.; });
  }

  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  std::unique_ptr<ScalarFunctionType> sigma_s() const override
  {
    return std::make_unique<GenericScalarFunctionType>(0, [](const DomainType& x, const XT::Common::Parameter&) {
      if (x[0] > 2)
        return 10.;
      else if (x[0] > 1)
        return 2.;
      else
        return 0.;
    });
  }

  // Q = 0.5 if 1 <= x <= 1.5, 0 else
  std::unique_ptr<ScalarFunctionType> Q() const override
  {
    return std::make_unique<GenericScalarFunctionType>(
        0, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] < 1 || x[0] > 1.5 ? 0. : 0.5; });
  }

protected:
  struct helper_base
  {
    // returns the numerator g of the left boundary value (see create_boundary_values)
    static RangeFieldType g(const RangeFieldType v)
    {
      return std::exp(-1e5 * (v - 1) * (v - 1));
    }

    static RangeFieldType dg(const RangeFieldType v)
    {
      return -2e5 * (v - 1) * g(v);
    }

    // returns the denominator <g> of the left boundary value (see create_boundary_values)
    static RangeFieldType denominator()
    {
      static RangeFieldType ret = 1 / 200. * std::sqrt(M_PI / 10) * std::erf(200 * std::sqrt(10));
      return ret;
    }

    // calculates integral from v_l to v_u of numerator g
    static RangeFieldType integral_1(RangeFieldType v_l, RangeFieldType v_u)
    {
      return 1 / 200. * std::sqrt(M_PI / 10)
             * (std::erf(100 * std::sqrt(10) * (v_u - 1)) - std::erf(100 * std::sqrt(10) * (v_l - 1)));
    }

    // calculates integral from v_l to v_u of v*g
    static RangeFieldType integral_2(RangeFieldType v_l, RangeFieldType v_u)
    {
      return integral_1(v_l, v_u) - 1. / 2e5 * (g(v_u) - g(v_l));
    }

    // calculates integral from v_l to v_u of v^2*g
    static RangeFieldType integral_3(RangeFieldType v_l, RangeFieldType v_u)
    {
      return 0.25e-10 * (dg(v_u) - dg(v_l)) + 2. * integral_2(v_l, v_u) + (0.5e-5 - 1.) * integral_1(v_l, v_u);
    }
  };


  template <class B, class anything = void>
  struct helper : public helper_base
  {
    using helper_base::denominator;
    using helper_base::g;

    static void get_left_boundary_values(const MomentBasisImp& basis_functions,
                                         const RangeFieldType& psi_vac,
                                         const bool is_mn_model,
                                         DynamicRangeType& ret)
    {
      // For the MN-Models, we have to use the quadrature also used in the optimization problem to guarantee
      // realizability of the boundary_values.
      // For the PN-Models, we do not have these issues and just use a very fine quadrature (which is not a performance
      // problem as the integration is only done once).
      std::fill(ret.begin(), ret.end(), 0.);
      const auto& quadratures =
          is_mn_model ? basis_functions.quadratures() : MomentBasisImp::gauss_lobatto_quadratures(100, 31);
      for (size_t ii = 0; ii < quadratures.size(); ++ii) {
        const auto& quadrature = quadratures[ii];
        for (const auto& quad_point : quadrature) {
          const auto& v = quad_point.position()[0];
          auto summand = basis_functions.evaluate(v, ii);
          summand *= g(v) * quad_point.weight();
          ret += summand;
        }
      }
      ret /= denominator();
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
    }

    static DynamicRangeType get_kinetic_boundary_flux(const MomentBasisImp& /*basis_functions*/,
                                                      const RangeFieldType& /*psi_vac*/,
                                                      const bool /*is_mn_model*/,
                                                      const DomainType& x,
                                                      const RangeFieldType& n,
                                                      const size_t dd,
                                                      const ThisType& problem)
    {
      return problem.kinetic_boundary_flux_from_quadrature(x, n, dd);
    }
  };

  template <class anything, EntropyType entropy>
  struct helper<HatFunctionMomentBasis<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, 1, entropy>, anything>
    : public helper_base
  {
    using helper_base::denominator;
    using helper_base::g;
    using helper_base::integral_1;
    using helper_base::integral_2;
    using helper_base::integral_3;

    static void get_left_boundary_values(const MomentBasisImp& basis_functions,
                                         const RangeFieldType psi_vac,
                                         const bool /*is_mn_model*/,
                                         DynamicRangeType& ret)
    {
      std::fill(ret.begin(), ret.end(), 0.);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        const auto& partitioning = basis_functions.partitioning();
        const auto mu_n = partitioning[nn];
        if (nn < dimRange - 1) {
          const auto mu_np1 = partitioning[nn + 1];
          ret[nn] +=
              1. / ((mu_n - mu_np1) * denominator()) * (integral_2(mu_n, mu_np1) - mu_np1 * integral_1(mu_n, mu_np1));
        }
        if (nn > 0) {
          const auto mu_nm1 = partitioning[nn - 1];
          ret[nn] +=
              1. / ((mu_n - mu_nm1) * denominator()) * (integral_2(mu_nm1, mu_n) - mu_nm1 * integral_1(mu_nm1, mu_n));
        }
      }
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
    } // ... get_left_boundary_values(...)

    static DynamicRangeType get_kinetic_boundary_flux(const MomentBasisImp& basis_functions,
                                                      const RangeFieldType& psi_vac,
                                                      const bool is_mn_model,
                                                      const DomainType& x,
                                                      const RangeFieldType& n,
                                                      const size_t dd,
                                                      const ThisType& problem)
    {
      if (!is_mn_model)
        DUNE_THROW(Dune::NotImplemented, "Only implemented for mn");
      if (x < 1.5) {
        DynamicRangeType ret(dimRange, 0.);
        const auto& partitioning = basis_functions.partitioning();
        for (size_t nn = 0; nn < dimRange; ++nn) {
          const auto& mu_n = partitioning[nn];
          if (nn < dimRange - 1) {
            const auto& mu_np1 = partitioning[nn + 1];
            if (mu_np1 > 0.) {
              const auto left_limit = mu_n > 0. ? mu_n : 0.;
              ret[nn] +=
                  1. / ((mu_n - mu_np1) * denominator())
                      * (integral_3(left_limit, mu_np1) - mu_np1 * integral_2(left_limit, mu_np1))
                  + psi_vac / 6. * std::pow(mu_np1 - left_limit, 2) * (mu_np1 + 2. * left_limit) / (mu_np1 - mu_n);
            } // if (mu_np1 > 0.)
          } // if (nn < dimRange -1)
          if (mu_n > 0.) {
            if (nn > 0) {
              const auto& mu_nm1 = partitioning[nn - 1];
              const auto left_limit = mu_nm1 > 0. ? mu_nm1 : 0.;
              ret[nn] += 1. / ((mu_n - mu_nm1) * denominator())
                             * (integral_3(left_limit, mu_n) - mu_nm1 * integral_2(left_limit, mu_n))
                         + psi_vac / 6.
                               * (3 * std::pow(mu_n, 2) * mu_nm1 - 3 * mu_nm1 * std::pow(left_limit, 2)
                                  - 2 * std::pow(mu_n, 3) + 2 * std::pow(left_limit, 3))
                               / (mu_nm1 - mu_n);
            } // if (nn > 0)
          } // if (mu_n > 0.)
        } // nn
        // unit outer normal
        assert(XT::Common::FloatCmp::eq(n, -1.));
        return ret;
      } else {
        return problem.kinetic_boundary_flux_from_quadrature(x, n, dd);
      }
    } // ... get_kinetic_boundary_flux(...)
  };

  template <class anything, EntropyType entropy>
  struct helper<PartialMomentBasis<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, 1, 1, entropy>, anything>
    : public helper_base
  {
    using helper_base::denominator;
    using helper_base::integral_1;
    using helper_base::integral_2;
    using helper_base::integral_3;

    static void get_left_boundary_values(const MomentBasisImp& basis_functions,
                                         const RangeFieldType psi_vac,
                                         const bool /*is_mn_model*/,
                                         DynamicRangeType& ret)
    {
      const auto& partitioning = basis_functions.partitioning();
      std::fill(ret.begin(), ret.end(), 0.);
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        ret[2 * ii] = integral_1(partitioning[ii], partitioning[ii + 1]) / denominator();
        ret[2 * ii + 1] = integral_2(partitioning[ii], partitioning[ii + 1]) / denominator();
      }
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
    }

    static DynamicRangeType get_kinetic_boundary_flux(const MomentBasisImp& basis_functions,
                                                      const RangeFieldType& psi_vac,
                                                      const bool is_mn_model,
                                                      const DomainType& x,
                                                      const RangeFieldType& n,
                                                      const size_t dd,
                                                      const ThisType& problem)
    {
      if (!is_mn_model)
        DUNE_THROW(Dune::NotImplemented, "Only implemented for mn");
      if (x < 1.5) {
        const auto& partitioning = basis_functions.partitioning();
        DynamicRangeType ret(dimRange, 0.);
        for (size_t ii = 0; ii < dimRange / 2; ++ii) {
          const auto& mu_i = partitioning[ii];
          const auto& mu_ip1 = partitioning[ii + 1];
          if (mu_ip1 > 0.) {
            const auto left_limit = mu_i > 0. ? mu_i : 0.;
            ret[2 * ii] = integral_2(left_limit, mu_ip1) / denominator()
                          + 0.5 * psi_vac * (std::pow(mu_ip1, 2) - std::pow(left_limit, 2));
            ret[2 * ii + 1] = integral_3(left_limit, partitioning[ii + 1]) / denominator()
                              + 1. / 3. * psi_vac * (std::pow(mu_ip1, 3) - std::pow(left_limit, 3));
          }
        } // ii
        assert(XT::Common::FloatCmp::eq(n, -1.));
        return ret;
      } else {
        return problem.kinetic_boundary_flux_from_quadrature(x, n, dd);
      }
    } // ... get_kinetic_boundary_flux(...)
  };

  using BaseType::basis_functions_;
  using BaseType::psi_vac_;
  const bool is_mn_model_;
}; // class SourceBeamPn<...>

template <class GV, class MomentBasis>
class SourceBeamMn : public SourceBeamPn<XT::Grid::extract_entity_t<GV>, MomentBasis>
{
  using BaseType = SourceBeamPn<XT::Grid::extract_entity_t<GV>, MomentBasis>;
  using ThisType = SourceBeamMn;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_grid_cfg;

  SourceBeamMn(const MomentBasis& basis_functions,
               const GV& grid_view,
               const RangeFieldType& psi_vac = 5e-7,
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const bool disable_realizability_check = false,
               const RangeFieldType tau = 1e-9)
    : BaseType(basis_functions, psi_vac, grid_cfg, true)
    , grid_view_(grid_view)
    , disable_realizability_check_(disable_realizability_check)
    , tau_(tau)
  {}

  static std::string static_id()
  {
    return "sourcebeammn";
  }

  std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_, tau_, disable_realizability_check_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
  const bool disable_realizability_check_;
  const RangeFieldType tau_;
}; // class SourceBeamMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
