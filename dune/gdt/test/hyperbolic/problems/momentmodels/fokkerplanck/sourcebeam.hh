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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_SOURCEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_SOURCEBEAM_HH

#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/legendre.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace FokkerPlanck {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class SourceBeamPn : public FokkerPlanckEquation<BasisfunctionImp, GridLayerImp, U_>
{
  typedef FokkerPlanckEquation<BasisfunctionImp, GridLayerImp, U_> BaseType;

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

  using BaseType::default_boundary_cfg;

  SourceBeamPn(const BasisfunctionType& basis_functions,
               const GridLayerType& grid_layer,
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, 6, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "sourcebeamfokkerplanckpn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[3.0]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{1, 1, 1, 1, 0, 0}),
                                  std::make_pair("T", std::vector<double>{0, 0, 2, 2, 10, 10}),
                                  std::make_pair("Q", std::vector<double>{0, 0, 1, 0, 0, 0}),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{4.0})});
  }

  // Boundary value of kinetic equation is delta(v-1) at x = 0 and psi_vac at x = 3,
  // so n-th component of boundary value has to be 0.5*\phi_n(1) at x = 0 and basis_integrated_n*psi_vac at x = 3.
  // Model with linear interpolating function 0.5*phi*(1-x/3) + basis_integrated*psi_vac * x/3
  virtual BoundaryValueType* create_boundary_values() const override
  {
    const auto basis_evaluated_at_one = basis_functions_.evaluate(DomainType(1));
    const auto basis_integrated = basis_functions_.integrated();
    return new ActualBoundaryValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          RangeType ret = basis_integrated;
          ret *= x[0] / 3. * psi_vac_;
          RangeType summand2 = basis_evaluated_at_one;
          summand2 *= (1 - x[0] / 3.) * 0.5;
          ret += summand2;
          return ret;
        },
        1);
  } // ... create_boundary_values()

protected:
  using BaseType::basis_functions_;
  using BaseType::psi_vac_;
}; // class SourceBeamPn<...>


} // namespace FokkerPlanck
} // namespace Problems


template <class G, class R = double>
class SourceBeamTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               Problems::KineticEquation<Problems::FokkerPlanck::
                                                             SourceBeamPn<Hyperbolic::Problems::
                                                                              LegendrePolynomials<double, double, 5>,
                                                                          typename G::LevelGridView,
                                                                          typename internal::
                                                                              DiscreteFunctionProvider<G,
                                                                                                       GDT::SpaceType::
                                                                                                           product_fv,
                                                                                                       0,
                                                                                                       R,
                                                                                                       6,
                                                                                                       1,
                                                                                                       GDT::Backends::
                                                                                                           gdt>::type>>>
{
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  typedef typename Hyperbolic::Problems::LegendrePolynomials<double, double, 5> BasisfunctionType;
  static const size_t dimRange = 6;
  static const size_t dimRangeCols = 1;
  typedef
      typename internal::DiscreteFunctionProvider<G, GDT::SpaceType::product_fv, 0, R, 6, 1, GDT::Backends::gdt>::type
          U;
  typedef typename Problems::
      KineticEquation<Problems::FokkerPlanck::SourceBeamPn<BasisfunctionType, typename G::LevelGridView, U>>
          ProblemType;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  SourceBeamTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
    , basis_functions_()
    , level_grid_view_(BaseType::level_view(0))
    , problem_(basis_functions_, level_grid_view_)
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
    out << "+======================================================================================================+"
           "\n"
        << "|+====================================================================================================+|"
           "\n"
        << "||  Testcase: Fokker-Planck SourceBeam                                                                ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|"
           "\n"
        << "||  domain = [0, 3]                                                                                   ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end()) + "]                                ||\n"
        << "||  flux = see http://dx.doi.org/10.1137/130934210 Section 6.5                                        ||\n"
        << "||  rhs = http://dx.doi.org/10.1137/130934210 Section 6.5                                             ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|"
           "\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const BasisfunctionType basis_functions_;
  const typename G::LevelGridView level_grid_view_;
  const ProblemType problem_;
}; // class SourceBeamTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_SOURCEBEAM_HH
