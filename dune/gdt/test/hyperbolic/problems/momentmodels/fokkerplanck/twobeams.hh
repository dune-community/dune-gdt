// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include "fokkerplanckequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace FokkerPlanck {


template <class BasisfunctionImp,
          class GridLayerImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class TwoBeamsPn : public FokkerPlanckEquation<BasisfunctionImp,
                                               GridLayerImp,
                                               EntityImp,
                                               DomainFieldImp,
                                               dimDomain,
                                               U_,
                                               RangeFieldImp,
                                               dimRange>
{
  typedef FokkerPlanckEquation<BasisfunctionImp,
                               GridLayerImp,
                               EntityImp,
                               DomainFieldImp,
                               dimDomain,
                               U_,
                               RangeFieldImp,
                               dimRange>
      BaseType;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  TwoBeamsPn(const BasisfunctionType& basis_functions,
             const GridLayerType& grid_layer,
             const QuadratureType& quadrature = default_quadrature(),
             const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
             const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, 1, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "twobeamspn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-0.5]";
    grid_config["upper_right"] = "[0.5]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1]";
    grid_config["num_quad_cells"] = "[20]";
    grid_config["quad_order"] = "50";
    return grid_config;
  }

  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{4}),
                                  std::make_pair("T", std::vector<double>{0}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.5}),
                                  std::make_pair("t_end", std::vector<double>{4.0}),
                                  std::make_pair("num_segments", std::vector<double>{1.})});
  }


  // boundary value of kinetic equation is 100*delta(v-1) at x = -0.5 and 100*delta(v+1) at x = 0.5,
  // so k-th component of boundary value has to be 50*\phi_k(1) at x = -0.5 and 50*\phi_k(-1) at x = 0.5.
  // Model with function(x) = 50*\phi_k(-1)*(x+0.5) + 50*\phi_k(1)*(0.5-x).
  virtual BoundaryValueType* create_boundary_values() const override
  {
    const auto basis_evaluated_at_one = basis_functions_.evaluate(DomainType(1));
    const auto basis_evaluated_at_minus_one = basis_functions_.evaluate(DomainType(-1));
    return new ActualBoundaryValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          RangeType ret = basis_evaluated_at_minus_one;
          ret *= (x[0] + 0.5) * 50;
          RangeType summand2 = basis_evaluated_at_one;
          summand2 *= (0.5 - x[0]) * 50;
          ret += summand2;
          return ret;
        },
        1);
  } // ... create_boundary_values()

protected:
  using BaseType::basis_functions_;
}; // class TwoBeamsPn<...>

template <class BasisfunctionType,
          class GridLayerType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class TwoBeamsMn : public TwoBeamsPn<BasisfunctionType,
                                     GridLayerType,
                                     EntityType,
                                     DomainFieldType,
                                     dimDomain,
                                     U_,
                                     RangeFieldType,
                                     dimRange>
{
  typedef TwoBeamsPn<BasisfunctionType,
                     GridLayerType,
                     EntityType,
                     DomainFieldType,
                     dimDomain,
                     U_,
                     RangeFieldType,
                     dimRange>
      BaseType;
  typedef TwoBeamsMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  using typename BaseType::QuadratureType;
  typedef EntropyBasedLocalFlux<BasisfunctionType,
                                GridLayerType,
                                EntityType,
                                DomainFieldType,
                                dimDomain,
                                U_,
                                RangeFieldType,
                                dimRange>
      ActualFluxType;


  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  TwoBeamsMn(const BasisfunctionType& basis_functions,
             const GridLayerType& grid_layer,
             const QuadratureType& quadrature = default_quadrature(),
             const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
             const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "twobeamsmn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_, quadrature_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
  using BaseType::quadrature_;
}; // class TwoBeamsMn<...>


} // namespace FokkerPlanck
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
