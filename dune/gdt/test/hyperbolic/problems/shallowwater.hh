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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>
#include <dune/xt/functions/lambda/global-function.hh>

#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

template <class E, class D, class U, size_t d = E::dimension>
class ShallowWater : public ProblemBase<E, D, d, U, typename U::RangeFieldType, 3>
{
  static_assert(d == 2, "This is the specialization for two dimension!");
  using ThisType = ShallowWater<E, D, U, d>;
  using BaseType = ProblemBase<E, D, d, U, typename U::RangeFieldType, 3>;

public:
  static const bool linear = false;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::ActualDirichletBoundaryValueType;
  using typename BaseType::ActualFluxType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::RhsType;
  using typename BaseType::StateRangeType;
  using ActualRhsType = typename XT::Functions::GlobalLambdaFluxFunction<U, 0, RangeFieldType, dimRange, 1>;

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0 0]";
    grid_config["upper_right"] = "[10 10]";
    grid_config["num_elements"] = "[100 100]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  using BaseType::default_boundary_cfg;

  ShallowWater(const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(create_flux(),
               create_rhs(),
               create_initial_values(grid_cfg),
               create_boundary_values(),
               grid_cfg,
               boundary_cfg,
               0.4,
               0.3,
               false)
  {}

  static std::string static_id()
  {
    return "ShallowWater";
  }

  static FluxType* create_flux()
  {
    return new ActualFluxType(
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          typename FluxType::RangeType ret{{u[1], u[2]},
                                           {std::pow(u[1], 2) / u[0] + 0.5 * std::pow(u[0], 2), u[1] * u[2] / u[0]},
                                           {u[1] * u[2] / u[0], std::pow(u[2], 2) / u[0] + 0.5 * std::pow(u[0], 2)}};
          return ret;
        },
        {},
        "shallow water flux",
        [](const XT::Common::Parameter&) { return 10; },
        [](const DomainType&, const StateRangeType&, const XT::Common::Parameter&) {
          return typename ActualFluxType::PartialXRangeType(0);
        },
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          typename ActualFluxType::PartialURangeType ret;
          ret[0] = {{0, 1, 0},
                    {-std::pow(u[1] / u[0], 2) + u[0], 2. * u[1] / u[0], 0.},
                    {-u[1] * u[2] / std::pow(u[0], 2), u[2] / u[0], u[1] / u[0]}};
          ret[1] = {{0, 0, 1},
                    {-u[1] * u[2] / std::pow(u[0], 2), u[2] / u[0], u[1] / u[0]},
                    {-std::pow(u[2] / u[0], 2) + u[0], 0, 2. * u[2] / u[0]}};
          return ret;
        });
  }

  static RhsType* create_rhs()
  {
    return new ActualRhsType(
        [](const DomainType& x, const StateRangeType& u, const XT::Common::Parameter& param) {
          if (x[0] > 6 && x[0] < 7 && x[1] > 6 && x[1] < 7 && param.get("t")[0] < 1)
            return RangeType{u[0], 0, 0};
          else
            return RangeType{0, 0, 0};
        },
        {},
        "shallow water rhs",
        [](const XT::Common::Parameter&) { return 10; });
  } // ... create_rhs(...)

  static InitialValueType* create_initial_values(const XT::Common::Configuration& grid_cfg)
  {
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg["upper_right"]);
    return new ActualInitialValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          return RangeType{std::exp(1 - std::pow(M_PI * (x[0] - upper_right[0] / 2), 2)
                                    - std::pow(M_PI * (x[1] - upper_right[1] / 2), 2))
                               + 1.,
                           0.,
                           0.};
        },
        10);
  } // ... create_initial_values()

  virtual BoundaryValueType* create_boundary_values()
  {
    return new ActualBoundaryValueType(std::make_unique<ActualDirichletBoundaryValueType>(
        [=](const DomainType&, const XT::Common::Parameter&) {
          return RangeType{1., 0., 0.};
        },
        0));
  } // ... create_boundary_values()
}; // class ShallowWater<...>


template <class E, class D, class U>
class ShallowWater<E, D, U, 1> : public ProblemBase<E, D, 1, U, typename U::RangeFieldType, 2>
{
  using ThisType = ShallowWater<E, D, U, 1>;
  using BaseType = ProblemBase<E, D, 1, U, typename U::RangeFieldType, 2>;

public:
  static const bool linear = false;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::ActualDirichletBoundaryValueType;
  using typename BaseType::ActualFluxType;
  using typename BaseType::ActualRhsType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using ActualInitialValueType = XT::Functions::
      CheckerboardFunction<E, D, dimDomain, RangeFieldType, dimRange, 1, ActualDirichletBoundaryValueType>;

  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;

  using typename BaseType::BoundaryValueType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using typename BaseType::RhsType;

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[10.0]";
    grid_config["num_elements"] = "[10]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  using BaseType::default_boundary_cfg;

  ShallowWater(const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(create_flux(),
               create_rhs(),
               create_initial_values(grid_cfg),
               create_boundary_values(),
               grid_cfg,
               boundary_cfg,
               0.4,
               3.,
               false)
  {}

  static std::string static_id()
  {
    return "ShallowWater";
  }

  static FluxType* create_flux()
  {
    return new ActualFluxType(
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          typename FluxType::RangeType ret;
          ret[0] = u[1];
          ret[1] = std::pow(u[1], 2) / u[0] + 0.5 * std::pow(u[0], 2);
          return ret;
        },
        {},
        "shallow water flux",
        [](const XT::Common::Parameter&) { return 10; },
        [](const DomainType&, const StateRangeType&, const XT::Common::Parameter&) {
          return typename ActualFluxType::PartialXRangeType(0);
        },
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          typename ActualFluxType::PartialURangeType ret;
          ret[0][0] = 0;
          ret[0][1] = 1;
          ret[1][0] = -std::pow(u[1] / u[0], 2) + u[0];
          ret[1][1] = 2. * u[1] / u[0];
          return ret;
        });
  }

  static RhsType* create_rhs()
  {
    return new ActualRhsType(FieldVector<MatrixType, 1>(MatrixType(0.)));
  } // ... create_rhs(...)

  static InitialValueType* create_initial_values(const XT::Common::Configuration& grid_cfg)
  {
    using LambdaFunctionType = typename ActualInitialValueType::LocalizableFunctionType;
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg["upper_right"]);
    const size_t num_regions = 5;
    FieldVector<size_t, dimDomain> num_segments(5);

    std::vector<LambdaFunctionType> initial_vals(num_regions,
                                                 LambdaFunctionType(
                                                     [](const DomainType&, const XT::Common::Parameter&) {
                                                       return RangeType{1., 0.};
                                                     },
                                                     0));

    initial_vals[2] = LambdaFunctionType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          return RangeType{1
                               + std::pow(x[0] - 4, 2) * std::pow(x[0] - 6, 2)
                                     * std::exp(2 - std::pow(x[0] - 4, 2) - std::pow(x[0] - 6, 2)),
                           0.};
        },
        10);
    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

  virtual BoundaryValueType* create_boundary_values()
  {
    return new ActualBoundaryValueType(XT::Grid::make_alldirichlet_boundaryinfo<IntersectionType>(),
                                       std::make_unique<ActualDirichletBoundaryValueType>(
                                           [=](const DomainType&, const XT::Common::Parameter&) {
                                             return RangeType{1., 0.};
                                           },
                                           0));
  } // ... create_boundary_values()
}; // class ShallowWater<...>


} // namespace Problems


// Test case for shallow water equations, see LeVeque, Finite Volume Methods for Hyperbolic Problems, 2002, Example 13.1
template <class G, class R = double>
class ShallowWaterTestCase
  : public Dune::GDT::Test::InstationaryTestCase<
        G,
        Problems::ShallowWater<typename G::template Codim<0>::Entity,
                               typename G::ctype,
                               typename internal::DiscreteFunctionProvider<G,
                                                                           GDT::SpaceType::product_fv,
                                                                           0,
                                                                           R,
                                                                           2,
                                                                           1,
                                                                           GDT::Backends::gdt,
                                                                           XT::LA::default_backend,
                                                                           XT::Grid::Layers::leaf,
                                                                           true

                                                                           >::type>>
{
  using E = typename G::template Codim<0>::Entity;
  using D = typename G::ctype;
  static const size_t d = G::dimension;

public:
  static const size_t dimRange = 2;
  static const size_t dimRangeCols = 1;
  using U = typename internal::DiscreteFunctionProvider<G,
                                                        GDT::SpaceType::product_fv,
                                                        0,
                                                        R,
                                                        2,
                                                        1,
                                                        GDT::Backends::gdt,
                                                        XT::LA::default_backend,
                                                        XT::Grid::Layers::leaf,
                                                        true>::type;
  using ProblemType = typename Problems::ShallowWater<E, D, U>;

private:
  using BaseType = typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType>;

public:
  using typename BaseType::GridType;

  ShallowWaterTestCase(const size_t num_refs = 2, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
  {}

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
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: Shallow Water                                           ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        << "||  domain = [0, 10]                                                  ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                   ||\n"
        << "||  flux = [u[1] u[1]*u[1]/u[0]+0.5*u[0]*u[0]]                        ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: solution on finest grid                       ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class ShallowWaterTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH
