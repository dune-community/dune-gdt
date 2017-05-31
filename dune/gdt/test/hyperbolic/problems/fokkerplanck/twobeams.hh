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

#include <dune/gdt/local/fluxes/entropybased.hh>

#include "kineticequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class TwoBeamsFokkerPlanckPn : public KineticFokkerPlanckEquation<BasisfunctionImp,
                                                                  EntityImp,
                                                                  DomainFieldImp,
                                                                  dimDomain,
                                                                  U_,
                                                                  RangeFieldImp,
                                                                  dimRange>
{
  typedef KineticFokkerPlanckEquation<BasisfunctionImp,
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

  using BaseType::default_boundary_cfg;

  TwoBeamsFokkerPlanckPn(const BasisfunctionType& basis_functions,
                         const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                         const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "twobeamspn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[1.0]";
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
                                  std::make_pair("num_elements", std::vector<double>{1.})});
  }


  // boundary value of kinetic equation is 100*delta(v-1) at x = 0 and 100*delta(v+1) at x = 1,
  // so k-th component of boundary value has to be 50*\phi_k(1) at x = 0 and 50*\phi_k(-1) at x = 1.
  // Model with function(x) = 50*\phi_k(-1)*x + 50*\phi_k(1)*(1-x).
  virtual BoundaryValueType* create_boundary_values() const override
  {
    const auto basis_evaluated_at_one = basis_functions_.evaluate(DomainType(1));
    const auto basis_evaluated_at_minus_one = basis_functions_.evaluate(DomainType(-1));
    return new ActualBoundaryValueType(
        [=](const DomainType& x) {
          RangeType ret = basis_evaluated_at_minus_one;
          ret *= x[0] * 50;
          RangeType summand2 = basis_evaluated_at_one;
          summand2 *= (1 - x[0]) * 50;
          ret += summand2;
          return ret;
        },
        1);
  } // ... create_boundary_values()

protected:
  using BaseType::basis_functions_;
}; // class TwoBeamsFokkerPlanckPn<...>

template <class GridViewType,
          class BasisfunctionType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class TwoBeamsFokkerPlanckMn : public TwoBeamsFokkerPlanckPn<BasisfunctionType,
                                                             EntityType,
                                                             DomainFieldType,
                                                             dimDomain,
                                                             U_,
                                                             RangeFieldType,
                                                             dimRange>
{
  typedef TwoBeamsFokkerPlanckPn<BasisfunctionType,
                                 EntityType,
                                 DomainFieldType,
                                 dimDomain,
                                 U_,
                                 RangeFieldType,
                                 dimRange>
      BaseType;
  typedef TwoBeamsFokkerPlanckMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  using typename BaseType::QuadratureType;
  typedef EntropyBasedLocalFlux<GridViewType,
                                BasisfunctionType,
                                EntityType,
                                DomainFieldType,
                                dimDomain,
                                U_,
                                RangeFieldType,
                                dimRange>
      ActualFluxType;


  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  TwoBeamsFokkerPlanckMn(const BasisfunctionType& basis_functions,
                         const GridViewType& grid_view,
                         const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                         const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {
  }

  static std::string static_id()
  {
    return "twobeamsmn";
  }

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(grid_view_, quadrature_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GridViewType& grid_view_;
  using BaseType::quadrature_;
}; // class TwoBeamsFokkerPlanckMn<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
