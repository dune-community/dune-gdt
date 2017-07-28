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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_BURGERS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_BURGERS_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>
#include <dune/xt/functions/lambda/global-function.hh>

#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class E, class D, size_t d, class U, class R, size_t r>
class Burgers : public ProblemBase<E, D, d, U, R, r>
{
  typedef Burgers<E, D, d, U, R, r> ThisType;
  typedef ProblemBase<E, D, d, U, R, r> BaseType;

public:
  static const bool linear = false;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  typedef typename XT::Functions::GlobalLambdaFluxFunction<U, 0, R, r, d> ActualFluxType;
  typedef typename XT::Functions::AffineFluxFunction<E, D, d, U, R, r, 1> ActualRhsType;
  typedef XT::Functions::GlobalLambdaFunction<E, D, d, R, r, 1> ActualBoundaryValueType;
  typedef ActualBoundaryValueType ActualInitialValueType;

  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  Burgers(const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
          const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(create_flux(),
               create_rhs(),
               create_initial_values(grid_cfg),
               create_boundary_values(),
               grid_cfg,
               boundary_cfg,
               dimDomain == 1 ? 0.5 : (dimDomain == 2 ? 0.1 : 0.05),
               1.,
               false)
  {
  }

  static std::string static_id()
  {
    return "Burgers";
  }

  static FluxType* create_flux()
  {
    return new ActualFluxType(
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          return typename FluxType::RangeType(0.5 * u[0] * u[0]);
        },
        {},
        "burgers_flux",
        [](const XT::Common::Parameter&) { return 2; },
        FieldVector<typename ActualFluxType::ColPartialXLambdaType, dimDomain>(
            [](const DomainType&, const StateRangeType&, const XT::Common::Parameter&) {
              return typename ActualFluxType::ColPartialXRangeType(0);
            }),
        FieldVector<typename ActualFluxType::ColPartialULambdaType, dimDomain>(
            [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
              typename ActualFluxType::ColPartialURangeType ret(0);
              for (size_t rr = 0; rr < dimRange; ++rr)
                ret[rr][0] = u[0];
              return ret;
            }));
  }

  static RhsType* create_rhs()
  {
    return new ActualRhsType(FieldVector<MatrixType, 1>(MatrixType(0.)));
  } // ... create_rhs(...)

  static InitialValueType* create_initial_values(const XT::Common::Configuration& /*grid_cfg*/)
  {
    return new ActualInitialValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          if (dimDomain == 1)
            return RangeType(std::sin(M_PI * x[0]));
          else
            return RangeType(1.0 / 40.0 * std::exp(1 - (2 * M_PI * x[0] - M_PI) * (2 * M_PI * x[0] - M_PI)
                                                   - (2 * M_PI * x[1] - M_PI)
                                                         * (2 * M_PI * x[1] - M_PI))); // bump, only in 2D or higher
        },
        10);
  } // ... create_initial_values()

  virtual BoundaryValueType* create_boundary_values()
  {
    return new ActualBoundaryValueType([=](const DomainType&, const XT::Common::Parameter&) { return 0; }, 0);
  } // ... create_boundary_values()
}; // class Burgers<...>


} // namespace Problems


template <class G, class R = double, size_t r = 1>
class BurgersTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               Problems::Burgers<typename G::template Codim<0>::Entity,
                                                 typename G::ctype,
                                                 G::dimension,
                                                 typename GDT::DiscreteFunctionProvider<G,
                                                                                        GDT::SpaceType::product_fv,
                                                                                        0,
                                                                                        R,
                                                                                        r,
                                                                                        1,
                                                                                        GDT::Backends::gdt>::type,
                                                 R,
                                                 r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  static const size_t dimRange = r;
  static const size_t dimRangeCols = 1;
  typedef typename GDT::DiscreteFunctionProvider<G, GDT::SpaceType::product_fv, 0, R, r, 1, GDT::Backends::gdt>::type U;
  typedef typename Problems::Burgers<E, D, d, U, R, r> ProblemType;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  BurgersTestCase(const size_t num_refs = (d == 1 ? 4 : 1), const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
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

  virtual std::bitset<d> periodic_directions() const override final
  {
    std::bitset<d> periodic_dirs;
    periodic_dirs.set();
    return periodic_dirs;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    const std::string domainstring = (d == 1)
                                         ? "||  domain = [0, 1]                                                   ||\n"
                                         : "||  domain = [0, 1] x [0, 1]                                          ||\n";
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: Burgers                                                 ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        << domainstring
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                     ||\n"
        << "||  flux = 0.5*u[0]^2                                                 ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: discrete solution on finest grid              ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class BurgersTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_BURGERS_HH
