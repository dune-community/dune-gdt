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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/pdelab/common/crossproduct.hh>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/math.hh>

#include "fokkerplanckequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace Fokkerplanck {


template <class BasisfunctionImp,
          class GridLayerImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class SourceBeamPn : public FokkerPlanckEquation<BasisfunctionImp,
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

  SourceBeamPn(const BasisfunctionType& basis_functions,
               const GridLayerType& grid_layer,
               const QuadratureType& quadrature = default_quadrature(),
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, 6, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "sourcebeamfokkerplanckpn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[3.0]";
    grid_config["num_elements"] = "[300]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{1, 1, 1, 1, 0, 0}),
                                  std::make_pair("sigma_s", std::vector<double>{0, 0, 2, 2, 10, 10}),
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
  using BaseType::quadrature_;
  using BaseType::psi_vac_;
}; // class SourceBeamPn<...>


} // namespace FokkerPlanck
} // namespace Problems


// template <class G,
//          class U,
//          class R = double,
//          size_t momentOrder = 5,
//          class B = Hyperbolic::Problems::LegendrePolynomials<typename G::ctype, typename G::ctype, momentOrder>>
// class SourceBeamTestCase
//    : public Dune::GDT::Test::
//          NonStationaryTestCase<G,
//                                Problems::Fokkerplanck::SourceBeamPn<B,
//    typename G::LevelGridView,
//                                                                   typename G::template Codim<0>::Entity,
//                                                                   typename G::ctype,
//                                                                   G::dimension,
//                                                                   U,
//                                                                   R,
//                                                                   momentOrder + 1>>
//{
//  typedef typename G::LevelGridView GV;
//  typedef typename G::template Codim<0>::Entity E;
//  typedef typename G::ctype D;

// public:
//  static const size_t d = G::dimension;
//  static_assert(d == 1, "Only implemented for dimension 1.");
//  typedef typename Problems::Fokkerplanck::SourceBeamPn<B, GV, E, D, d, U, R, momentOrder + 1> ProblemType;
//  static const size_t dimRange = ProblemType::dimRange;
//  static const size_t dimRangeCols = 1;

// private:
//  typedef typename Dune::GDT::Test::NonStationaryTestCase<G, ProblemType> BaseType;

// public:
//  using typename BaseType::GridType;
//  using typename BaseType::SolutionType;

//  SourceBeamTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
//    : BaseType(
//          divide_t_end_by, XT::Grid::make_cube_grid<GridType>(ProblemType::default_grid_config()).grid_ptr(),
//          num_refs)
//    , problem_(B())
//  {
//  }

//  virtual const ProblemType& problem() const override final
//  {
//    return problem_;
//  }

//  virtual bool provides_exact_solution() const override final
//  {
//    return false;
//  }

//  virtual void print_header(std::ostream& out = std::cout) const override final
//  {
//    out <<
//    "+======================================================================================================+\n"
//        <<
//        "|+====================================================================================================+|\n"
//        << "||  Testcase: Fokker-Planck SourceBeam ||\n"
//        <<
//        "|+----------------------------------------------------------------------------------------------------+|\n"
//        << "||  domain = [0, 3] ||\n"
//        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end()) + "] ||\n"
//        << "||  flux = see http://dx.doi.org/10.1137/130934210 Section 6.5 ||\n"
//        << "||  rhs = http://dx.doi.org/10.1137/130934210 Section 6.5 ||\n"
//        << "||  reference solution: discrete solution on finest grid ||\n"
//        <<
//        "|+====================================================================================================+|\n"
//        << "+======================================================================================================+"
//        << std::endl;
//  }

// private:
//  const ProblemType problem_;
//}; // class SourceBeamTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
