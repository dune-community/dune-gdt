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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORT_CHECKERBOARD_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORT_CHECKERBOARD_HH

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


/**
 * In 2D, this is the Testcase for the Boltzmann equation in two dimensions,
 * see Section 4.2 of Brunner, Holloway, "Two-dimensional time dependent Riemann solvers for neutron transport", Journal
 * of Computational Physics, Volume 210, Issue 1, 2005
 * http://dx.doi.org/10.1016/j.jcp.2005.04.011
 * The 3D version is a straightforward generalization of the setup to 3 dimensions.
 * */
template <class BasisfunctionImp, class GridLayerImp, class U_>
class CheckerboardPn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>
{
  typedef KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_> BaseType;

public:
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using BaseType::dimDomain;

  static XT::Common::Configuration default_boundary_cfg()
  {
    XT::Common::Configuration boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  using BaseType::default_quadrature;

  CheckerboardPn(const BasisfunctionType& basis_functions,
                 const GridLayerType& grid_layer,
                 const QuadratureType& quadrature = default_quadrature(),
                 const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                 const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, {7, 7, 7}, grid_cfg, boundary_cfg, 1e-8 / (4 * M_PI))
  {
  }

  static std::string static_id()
  {
    size_t domainDim = dimDomain; // avoid linker error for dimDomain
    return "checkerboard" + XT::Common::to_string(domainDim) + "d_pn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0 0.0 0.0]";
    grid_config["upper_right"] = "[7.0 7.0 7.0]";
    grid_config["num_elements"] = "[14 14 14]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  // see the Checkerboard test case in http://www.sciencedirect.com/science/article/pii/S0021999105002275?np=y
  // q = 0 except in the center where q = 1. Sigma_s = Sigma_t = 1 in scattering regions, Sigma_s = 0, Sigma_t
  // = 10 in absorbing regions. Center is also a scattering region.
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", create_sigma_a()),
                                  std::make_pair("sigma_s", create_sigma_s()),
                                  std::make_pair("Q", create_Q()),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{3.2})});
  }

protected:
  static std::vector<double> create_sigma_a()
  {
    std::vector<double> sigma_a(std::pow(7, dimDomain));
    if (dimDomain == 2)
      for (size_t row = 0; row < 7; ++row)
        for (size_t col = 0; col < 7; ++col)
          sigma_a[7 * row + col] = is_absorbing(row, col) ? 10. : 0.;
    else if (dimDomain == 3)
      for (size_t plane = 0; plane < 7; ++plane)
        for (size_t row = 0; row < 7; ++row)
          for (size_t col = 0; col < 7; ++col)
            sigma_a[49 * plane + 7 * row + col] = is_absorbing(plane, row, col) ? 10 : 0.;
    return sigma_a;
  }

  static std::vector<double> create_sigma_s()
  {
    std::vector<double> sigma_s(std::pow(7, dimDomain));
    if (dimDomain == 2)
      for (size_t row = 0; row < 7; ++row)
        for (size_t col = 0; col < 7; ++col)
          sigma_s[7 * row + col] = is_absorbing(row, col) ? 0 : 1.;
    else if (dimDomain == 3)
      for (size_t plane = 0; plane < 7; ++plane)
        for (size_t row = 0; row < 7; ++row)
          for (size_t col = 0; col < 7; ++col)
            sigma_s[49 * plane + 7 * row + col] = is_absorbing(plane, row, col) ? 0 : 1.;
    return sigma_s;
  }

  static std::vector<double> create_Q()
  {
    std::vector<double> Q(std::pow(7, dimDomain), 0.);
    if (dimDomain == 2)
      Q[3 * 7 + 3] = 1. / std::sqrt(4. * M_PI);
    else if (dimDomain == 3)
      Q[3 * 7 * 7 + 3 * 7 + 3] = 1. / std::sqrt(4. * M_PI);
    return Q;
  }

  static bool is_absorbing(const size_t row, const size_t col)
  {
    return (row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
           || ((row == 3 || row == 5) && (col == 1 || col == 5));
  }

  static bool is_absorbing(size_t plane, size_t row, size_t col)
  {
    DXT_ASSERT(plane < 7 && row < 7 && col < 7);
    if (plane == 0 || plane == 6)
      return false;
    if (plane == 2 || plane == 4)
      return !((row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
               || ((row == 3 || row == 5) && (col == 1 || col == 5)))
             && (row != 0 && row != 6) && (col != 0 && col != 6) && !(row == 5 && col == 3) && !(row == 3 && col == 3);
    if (plane == 3 || plane == 1 || plane == 5)
      return (row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
             || ((row == 3 || row == 5) && (col == 1 || col == 5))
             || (plane != 3 && col == 3 && (row == 3 || row == 5));
    return false;
  }
}; // class CheckerboardPn<...>

template <class BasisfunctionType, class GridLayerType, class U_>
class CheckerboardMn : public CheckerboardPn<BasisfunctionType, GridLayerType, U_>
{
  typedef CheckerboardPn<BasisfunctionType, GridLayerType, U_> BaseType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef GDT::EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, U_> ActualFluxType;
  using typename BaseType::QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  CheckerboardMn(const BasisfunctionType& basis_functions,
                 const GridLayerType& grid_layer,
                 const QuadratureType& quadrature,
                 const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                 const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "checkerboardmn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_, quadrature_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
  using BaseType::quadrature_;
}; // class CheckerboardMn<...>


} // namespace KineticTransport
} // namespace Problems


template <class G, class R = double, size_t momentOrder = 1>
class CheckerboardTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               Problems::
                                   KineticEquation<Problems::KineticTransport::
                                                       CheckerboardPn<Hyperbolic::Problems::
                                                                          SphericalHarmonics<double,
                                                                                             double,
                                                                                             momentOrder,
                                                                                             G::dimension,
                                                                                             true>,
                                                                      typename G::LevelGridView,
                                                                      typename internal::
                                                                          DiscreteFunctionProvider<G,
                                                                                                   GDT::SpaceType::
                                                                                                       product_fv,
                                                                                                   0,
                                                                                                   R,
                                                                                                   Hyperbolic::Problems::
                                                                                                       SphericalHarmonics<double,
                                                                                                                          double,
                                                                                                                          momentOrder,
                                                                                                                          G::dimension,
                                                                                                                          true>::
                                                                                                           dimRange,
                                                                                                   1,
                                                                                                   GDT::Backends::gdt>::
                                                                              type>>>
{
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef Hyperbolic::Problems::SphericalHarmonics<double, double, momentOrder, G::dimension, true> BasisfunctionType;

public:
  typedef Problems::
      KineticEquation<Problems::KineticTransport::
                          CheckerboardPn<BasisfunctionType,
                                         typename G::LevelGridView,
                                         typename internal::
                                             DiscreteFunctionProvider<G,
                                                                      GDT::SpaceType::product_fv,
                                                                      0,
                                                                      R,
                                                                      Hyperbolic::Problems::
                                                                          SphericalHarmonics<double,
                                                                                             double,
                                                                                             momentOrder,
                                                                                             G::dimension,
                                                                                             true>::dimRange,
                                                                      1,
                                                                      GDT::Backends::gdt>::type>>
          ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;

  CheckerboardTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
    , problem_(BasisfunctionType(), BaseType::level_view(0))
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return false;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+======================================================================================================+\n"
        << "|+====================================================================================================+|\n"
        << "||  Testcase: Boltzmann 2D Checkerboard                                                               ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [0, 7] x [0, 7]                                                                          ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                                                  ||\n"
        << "||  flux = see http://dx.doi.org/10.1016/j.jcp.2005.04.011 Section 4.1                                ||\n"
        << "||  rhs = see http://dx.doi.org/10.1016/j.jcp.2005.04.011 Section 4.1                                 ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const ProblemType problem_;
}; // class CheckerboardTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORT_CHECKERBOARD_HH
