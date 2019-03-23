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


template <class E, class BasisfunctionImp>
class SourceBeamPn : public KineticTransportEquationBase<E, BasisfunctionImp>
{
  using BaseType = KineticTransportEquationBase<E, BasisfunctionImp>;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::FluxType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::GenericScalarFunctionType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::ScalarFunctionType;

  using BaseType::default_boundary_cfg;

  SourceBeamPn(const BasisfunctionType& basis_functions,
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
  virtual std::unique_ptr<BoundaryValueType> boundary_values() const override final
  {
    return std::make_unique<GenericFunctionType>(1, [&](const DomainType& x, const XT::Common::Parameter&) {
      if (x[0] < 1.5) {
        static auto ret = helper<BasisfunctionType>::get_left_boundary_values(basis_functions_, psi_vac_, is_mn_model_);
        return ret;
      } else {
        auto ret = basis_functions_.integrated();
        ret *= psi_vac_;
        return ret;
      }
    });
  } // ... boundary_values()

  RangeReturnType left_boundary_value() const
  {
    return helper<BasisfunctionType>::get_left_boundary_values(basis_functions_, psi_vac_, is_mn_model_);
  }

  virtual RangeFieldType t_end() const override
  {
    return 2.5;
  }

  // sigma_a = 1 if x <= 2, 0 else
  virtual std::unique_ptr<ScalarFunctionType> sigma_a() const override
  {
    return std::make_unique<GenericScalarFunctionType>(
        0, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] > 2 ? 0. : 1.; });
  }

  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  virtual std::unique_ptr<ScalarFunctionType> sigma_s() const override
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
  virtual std::unique_ptr<ScalarFunctionType> Q() const override
  {
    return std::make_unique<GenericScalarFunctionType>(
        0, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] < 1 || x[0] > 1.5 ? 0. : 0.5; });
  }

protected:
  struct helper_base
  {
    // returns the numerator g of the left boundary value (see create_boundary_values)
    static RangeFieldType numerator(const RangeFieldType v)
    {
      return std::exp(-1e5 * (v - 1) * (v - 1));
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
      return integral_1(v_l, v_u) - 1. / 2e5 * (numerator(v_u) - numerator(v_l));
    }
  };


  template <class B, class anything = void>
  struct helper : public helper_base
  {
    using helper_base::denominator;
    using helper_base::numerator;

    static DynamicRangeType get_left_boundary_values(const BasisfunctionImp& basis_functions,
                                                     const RangeFieldType& psi_vac,
                                                     const bool is_mn_model)
    {
      DynamicRangeType ret(dimRange, 0);
      // For the MN-Models, we have to use the quadrature also used in the optimization problem to guarantee
      // realizability of the boundary_values.
      // For the PN-Models, we do not have these issues and just use a very fine quadrature (which is not a performance
      // problem as the integration is only done once).
      const auto& quadratures =
          is_mn_model ? basis_functions.quadratures() : BasisfunctionImp::gauss_lobatto_quadratures(100, 31);
      for (size_t ii = 0; ii < quadratures.size(); ++ii) {
        const auto& quadrature = quadratures[ii];
        for (const auto& quad_point : quadrature) {
          const auto& v = quad_point.position()[0];
          auto summand = basis_functions.evaluate(v, ii);
          summand *= numerator(v) * quad_point.weight();
          ret += summand;
        }
      }
      ret /= denominator();
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
      return ret;
    }
  };

  template <class anything>
  struct helper<HatFunctionMomentBasis<DomainFieldType, dimDomain, RangeFieldType, dimRange>, anything>
    : public helper_base
  {
    using helper_base::denominator;
    using helper_base::integral_1;
    using helper_base::numerator;

    static DynamicRangeType get_left_boundary_values(const BasisfunctionImp& basis_functions,
                                                     const RangeFieldType psi_vac,
                                                     const bool /*is_mn_model*/)
    {
      DynamicRangeType ret(dimRange, 0);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        const auto& triangulation = basis_functions.triangulation();
        const auto vn = triangulation[nn];
        if (nn < dimRange - 1) {
          const auto vnp = triangulation[nn + 1];
          ret[nn] += 1. / ((vn - vnp) * denominator())
                     * ((1 - vnp) * integral_1(vn, vnp) - 1. / 2e5 * (numerator(vnp) - numerator(vn)));
        }
        if (nn > 0) {
          const auto vnm = triangulation[nn - 1];
          ret[nn] += 1. / ((vn - vnm) * denominator())
                     * ((1 - vnm) * integral_1(vnm, vn) - 1. / 2e5 * (numerator(vn) - numerator(vnm)));
        }
      }
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
      return ret;
    }
  };

  template <class anything>
  struct helper<PartialMomentBasis<DomainFieldType, dimDomain, RangeFieldType, dimRange>, anything> : public helper_base
  {
    using helper_base::denominator;
    using helper_base::integral_1;
    using helper_base::integral_2;

    static DynamicRangeType get_left_boundary_values(const BasisfunctionImp& basis_functions,
                                                     const RangeFieldType psi_vac,
                                                     const bool /*is_mn_model*/)
    {
      const auto& triangulation = basis_functions.triangulation();
      DynamicRangeType ret(dimRange, 0);
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        ret[2 * ii] = integral_1(triangulation[ii], triangulation[ii + 1]) / denominator();
        ret[2 * ii + 1] = integral_2(triangulation[ii], triangulation[ii + 1]) / denominator();
      }
      // add small vacuum concentration to move away from realizable boundary
      ret += basis_functions.integrated() * psi_vac;
      return ret;
    }
  };

  using BaseType::basis_functions_;
  using BaseType::psi_vac_;
  const bool is_mn_model_;
}; // class SourceBeamPn<...>

template <class E, class BasisfunctionType>
class SourceBeamMn : public SourceBeamPn<E, BasisfunctionType>
{
  using BaseType = SourceBeamPn<E, BasisfunctionType>;
  using ThisType = SourceBeamMn;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeReturnType;
  using ActualFluxType = EntropyBasedLocalFlux<BasisfunctionType>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  SourceBeamMn(const BasisfunctionType& basis_functions,
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg, true)
  {}

  static std::string static_id()
  {
    return "sourcebeammn";
  }

  virtual std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
}; // class SourceBeamMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
