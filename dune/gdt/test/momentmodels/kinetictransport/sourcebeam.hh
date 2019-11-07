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

  using BaseType::default_boundary_cfg;

  SourceBeamPn(const MomentBasis& basis_functions,
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
               const bool is_mn_model = false)
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
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
    return std::make_unique<GenericFunctionType>(1, [&](const DomainType& x, const XT::Common::Parameter&) {
      if (x[0] < 1.5) {
        static auto ret = helper<MomentBasis>::get_left_boundary_values(basis_functions_, psi_vac_, is_mn_model_);
        return ret;
      } else {
        auto ret = basis_functions_.integrated();
        ret *= psi_vac_;
        return ret;
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
          return (val > this->psi_vac_) ? val : this->psi_vac_;
        };
    };
  }

  virtual RangeReturnType
  kinetic_boundary_flux(const DomainType& x, const RangeFieldType& n, const size_t dd) const override final
  {
    return helper<MomentBasis>::get_kinetic_boundary_flux(basis_functions_, psi_vac_, is_mn_model_, x, n, dd, *this);
  }

  RangeReturnType left_boundary_value() const
  {
    return helper<MomentBasis>::get_left_boundary_values(basis_functions_, psi_vac_, is_mn_model_);
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

    static DynamicRangeType get_left_boundary_values(const MomentBasisImp& basis_functions,
                                                     const RangeFieldType& psi_vac,
                                                     const bool is_mn_model)
    {
      DynamicRangeType ret(dimRange, 0);
      // For the MN-Models, we have to use the quadrature also used in the optimization problem to guarantee
      // realizability of the boundary_values.
      // For the PN-Models, we do not have these issues and just use a very fine quadrature (which is not a performance
      // problem as the integration is only done once).
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
      return ret;
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

    static DynamicRangeType get_left_boundary_values(const MomentBasisImp& basis_functions,
                                                     const RangeFieldType psi_vac,
                                                     const bool /*is_mn_model*/)
    {
      DynamicRangeType ret(dimRange, 0);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        const auto& partitioning = basis_functions.partitioning();
        const auto vn = partitioning[nn];
        if (nn < dimRange - 1) {
          const auto vnp = partitioning[nn + 1];
          ret[nn] += 1. / ((vn - vnp) * denominator()) * (integral_2(vn, vnp) - vnp * integral_1(vn, vnp));
        }
        if (nn > 0) {
          const auto vnm = partitioning[nn - 1];
          ret[nn] += 1. / ((vn - vnm) * denominator()) * (integral_2(vnm, vn) - vnm * integral_1(vnm, vn));
        }
      }
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
      return ret;
    } // ... get_left_boundary_values(...)

    static DynamicRangeType get_kinetic_boundary_flux(const MomentBasisImp& basis_functions,
                                                      const RangeFieldType& /*psi_vac*/,
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
          const auto vn = partitioning[nn];
          if (nn < dimRange - 1) {
            const auto vnp = partitioning[nn + 1];
            if (vnp > 0.) {
              const auto left_limit = vn > 0. ? vn : 0.;
              ret[nn] +=
                  1. / ((vn - vnp) * denominator()) * (integral_3(left_limit, vnp) - vnp * integral_2(left_limit, vnp));
            } // if (vnp > 0.)
          } // if (nn < dimRange -1)
          if (vn > 0.) {
            if (nn > 0) {
              const auto vnm = partitioning[nn - 1];
              const auto left_limit = vnm > 0. ? vnm : 0.;
              ret[nn] +=
                  1. / ((vn - vnm) * denominator()) * (integral_2(left_limit, vn) - vnm * integral_1(left_limit, vn));
            } // if (nn > 0)
          } // if (vn > 0.)
        } // nn
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

    static DynamicRangeType get_left_boundary_values(const MomentBasisImp& basis_functions,
                                                     const RangeFieldType psi_vac,
                                                     const bool /*is_mn_model*/)
    {
      const auto& partitioning = basis_functions.partitioning();
      DynamicRangeType ret(dimRange, 0);
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        ret[2 * ii] = integral_1(partitioning[ii], partitioning[ii + 1]) / denominator();
        ret[2 * ii + 1] = integral_2(partitioning[ii], partitioning[ii + 1]) / denominator();
      }
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
      return ret;
    }

    static DynamicRangeType get_kinetic_boundary_flux(const MomentBasisImp& basis_functions,
                                                      const RangeFieldType& /*psi_vac*/,
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
          if (partitioning[ii + 1] > 0.) {
            const auto left_limit = partitioning[ii] > 0. ? partitioning[ii] : 0.;
            ret[2 * ii] = integral_2(left_limit, partitioning[ii + 1]) / denominator();
            ret[2 * ii + 1] = integral_3(left_limit, partitioning[ii + 1]) / denominator();
          }
        } // ii
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
  using typename BaseType::RangeReturnType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  SourceBeamMn(const MomentBasis& basis_functions,
               const GV& grid_view,
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg, true)
    , grid_view_(grid_view)
  {}

  static std::string static_id()
  {
    return "sourcebeammn";
  }

  std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
}; // class SourceBeamMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
