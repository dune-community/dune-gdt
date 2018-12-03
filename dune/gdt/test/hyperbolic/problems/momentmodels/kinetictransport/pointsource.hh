// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class PointSourcePn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>
{
  using BaseType = KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>;

public:
  using typename BaseType::ActualDirichletBoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::GridLayerType;
  using typename BaseType::InitialValueType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;

  using BaseType::default_boundary_cfg;

  PointSourcePn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, {1, 1, 1}, grid_cfg, boundary_cfg, 1e-4 / (4 * M_PI))
  {}

  static std::string static_id()
  {
    return "pointsourcepn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-1 -1 -1]";
    grid_config["upper_right"] = "[1 1 1]";
    grid_config["num_elements"] = "[10 10 10]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  // sigma_a = 0, sigma_s = 1, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{1}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.49 * 1 / std::sqrt(3)}),
                                  std::make_pair("t_end", std::vector<double>{0.45})});
  }

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-||x||^2/(2*sigma^2)).
  // Thus the initial value for the moments is basis_integrated * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-||x||^2/(2*sigma^2))).
  virtual InitialValueType* create_initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    static const double sigma = 0.03;
    RangeType basis_integrated = basis_functions_.integrated();
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;

    initial_vals.emplace_back(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          static const auto first_factor = 1. / (4 * M_PI * std::pow(M_PI * sigma, 3));
          static const auto second_factor = 1. / (M_PI * std::pow(sigma, 2));
          return basis_integrated * std::max(first_factor * std::exp(-x.two_norm2() * second_factor), psi_vac_);
        },
        61);

    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_cfg_;
  using BaseType::num_segments_;
  using BaseType::psi_vac_;
}; // class PointSourcePn<...>

template <class BasisfunctionType, class GridLayerType, class U_>
class PointSourceMn : public PointSourcePn<BasisfunctionType, GridLayerType, U_>
{
  using BaseType = PointSourcePn<BasisfunctionType, GridLayerType, U_>;
  using ThisType = PointSourceMn;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  using ActualFluxType = GDT::EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, U_>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  PointSourceMn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, grid_cfg, boundary_cfg)
  {}

  static std::string static_id()
  {
    return "pointsourcemn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
}; // class PointSourceMn<...>


} // namespace KineticTransport
} // namespace Problems


template <class G,
          class R = double,
          size_t rangeDim = 6,
          class B = HatFunctionMomentBasis<typename G::ctype, 3, typename G::ctype, rangeDim, 1, 3>>
class PointSourceTestCase
  : public Dune::GDT::Test::InstationaryTestCase<
        G,
        typename Hyperbolic::Problems::KineticEquation<typename Problems::KineticTransport::PointSourcePn<
            B,
            typename G::LeafGridLayer,
            typename internal::DiscreteFunctionProvider<G, GDT::SpaceType::product_fv, 0, R, 6, 1, GDT::Backends::gdt>::
                type>>>
{
  using D = typename G::ctype;

public:
  static const size_t d = G::dimension;
  static_assert(d == 3, "Only implemented for dimension 3.");
  using ProblemType = typename Hyperbolic::Problems::KineticEquation<typename Problems::KineticTransport::PointSourcePn<
      B,
      typename G::LeafGridLayer,
      DiscreteFunction<FvProductSpace<typename G::LeafGridLayer, double, rangeDim, 1>,
                       typename Dune::XT::LA::Container<double, XT::LA::default_sparse_backend>::VectorType>>>;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  using BaseType = typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType>;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;

  PointSourceTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
    , problem_(B())
  {}

  //  virtual const ProblemType& problem() const override final
  //  {
  //    return problem_;
  //  }

  virtual bool provides_exact_solution() const override final
  {
    return false;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+======================================================================================================+\n"
        << "|+====================================================================================================+|\n"
        << "||  Testcase: PointSource Pn                                                                          ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [-0.5, 0.5]^3                                                                            ||\n"
        << "||  time = [0, 0.45]                                                                                  ||\n"
        << "||  flux = see http://dx.doi.org/10.1137/130934210 Section 6.5                                        ||\n"
        << "||  rhs = http://dx.doi.org/10.1137/130934210 Section 6.5                                             ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const ProblemType problem_;
}; // class PointSourceTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
