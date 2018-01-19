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

#include "kinetictransportequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp, class GridLayerImp, class U_, bool linear = true>
class PointSourcePn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_, linear>
{
  typedef KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_, linear> BaseType;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  PointSourcePn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const QuadratureType& quadrature = default_quadrature(),
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, {1, 1, 1}, grid_cfg, boundary_cfg)
  {
  }

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
    grid_config["num_elements"] = "[4 4 4]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  // sigma_a = 0, sigma_s = 1, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{1}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.3}),
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

    //    initial_vals.emplace_back(
    //        [=](const DomainType& x) {
    //          auto ret = basis_integrated;
    //          ret *= psi_vac_ + 1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm() / (2. * sigma *
    //          sigma));
    //          return ret;
    //        },
    //        50);

    initial_vals.emplace_back(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = basis_integrated;
          //          ret *= std::max(1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm2() / (2. * sigma *
          //          sigma)),
          //                          1e-4 / (4. * M_PI));
          ret *= std::max(1. / (4 * M_PI * std::pow(M_PI * sigma, 3))
                              * std::exp(-x.two_norm2() / (M_PI * std::pow(sigma, 2))),
                          1e-4 / (4. * M_PI));
          return ret;
        },
        50);

    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::grid_cfg_;
  using BaseType::basis_functions_;
  using BaseType::num_segments_;
  using BaseType::psi_vac_;
}; // class PointSourcePn<...>

template <class BasisfunctionType, class GridLayerType, class U_>
class PointSourceMn : public PointSourcePn<BasisfunctionType, GridLayerType, U_, false /*nonlinear*/>
{
  typedef PointSourcePn<BasisfunctionType, GridLayerType, U_, false> BaseType;
  typedef PointSourceMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef GDT::EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, U_> ActualFluxType;
  using typename BaseType::QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  PointSourceMn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const QuadratureType& quadrature,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "pointsourcemn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_, quadrature_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
  using BaseType::quadrature_;
}; // class PointSourceMn<...>


} // namespace KineticTransport
} // namespace Problems


template <class G,
          class R = double,
          size_t rangeDim = 6,
          class B = Hyperbolic::Problems::HatFunctions<typename G::ctype, 3, typename G::ctype, rangeDim, 1, 3>>
class PointSourceTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               typename Hyperbolic::Problems::KineticEquation<
                                   typename Problems::KineticTransport::
                                       PointSourcePn<B,
                                                     typename G::LeafGridLayer,
                                                     typename internal::
                                                         DiscreteFunctionProvider<G,
                                                                                  GDT::SpaceType::product_fv,
                                                                                  0,
                                                                                  R,
                                                                                  6,
                                                                                  1,
                                                                                  GDT::Backends::gdt>::type>>>
{
  typedef typename G::ctype D;

public:
  static const size_t d = G::dimension;
  static_assert(d == 3, "Only implemented for dimension 3.");
  typedef typename Hyperbolic::Problems::KineticEquation<
      typename Problems::KineticTransport::
          PointSourcePn<B,
                        typename G::LeafGridLayer,
                        DiscreteFunction<FvProductSpace<typename G::LeafGridLayer, double, rangeDim, 1>,
                                         typename Dune::XT::LA::Container<double,
                                                                          XT::LA::default_sparse_backend>::VectorType>>>
      ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;

  PointSourceTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
    , problem_(B())
  {
  }

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
