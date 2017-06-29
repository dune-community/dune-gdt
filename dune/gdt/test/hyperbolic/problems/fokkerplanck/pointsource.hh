// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/test/instationary-eocstudy.hh>

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
class PointSourcePn : public KineticTransportEquation<BasisfunctionImp,
                                                      EntityImp,
                                                      DomainFieldImp,
                                                      dimDomain,
                                                      U_,
                                                      RangeFieldImp,
                                                      dimRange>
{
  typedef KineticTransportEquation<BasisfunctionImp, EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange>
      BaseType;

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
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  PointSourcePn(const BasisfunctionType& basis_functions,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
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
    //    grid_config["lower_left"] = "[-0.5 -0.5 -0.5]";
    //    grid_config["upper_right"] = "[0.5 0.5 0.5]";
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
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{0.45}),
                                  std::make_pair("num_segments", std::vector<double>{1., 1., 1.})});
  }

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the moments is basis_integrated * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  virtual InitialValueType* create_initial_values() const
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const FieldVector<size_t, 3> num_segments = get_num_segments(parameters());
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
        [=](const DomainType& x) {
          auto ret = basis_integrated;
          ret *= std::max(1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm2() / (2. * sigma * sigma)),
                          1e-4 / (4. * M_PI));
          return ret;
        },
        50);

    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::get_num_segments;

  using BaseType::grid_cfg_;
  using BaseType::basis_functions_;
  using BaseType::psi_vac_;
}; // class PointSourcePn<...>

template <class GridViewType,
          class BasisfunctionType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class PointSourceMn
    : public PointSourcePn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange>
{
  typedef PointSourcePn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange>
      BaseType;
  typedef PointSourceMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef GDT::EntropyBasedLocalFlux<GridViewType,
                                     BasisfunctionType,
                                     EntityType,
                                     DomainFieldType,
                                     dimDomain,
                                     U_,
                                     RangeFieldType,
                                     dimRange>
      ActualFluxType;
  using typename BaseType::QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  PointSourceMn(const BasisfunctionType& basis_functions,
                const QuadratureType& quadrature,
                const GridViewType& grid_view,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
    , quadrature_(quadrature)
  {
  }

  static std::string static_id()
  {
    return "pointsourcemn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(grid_view_, quadrature_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GridViewType& grid_view_;
  const QuadratureType& quadrature_;
}; // class PointSourceMn<...>


} // namespace Problems


template <class G,
          class R = double,
          size_t rangeDim = 6,
          class B = Hyperbolic::Problems::HatFunctions<typename G::ctype, 3, typename G::ctype, rangeDim, 1, 3>>
class PointSourceTestCase
    : public Dune::GDT::Test::
          NonStationaryTestCase<G,
                                typename Hyperbolic::Problems::KineticEquation<
                                    typename Problems::
                                        PointSourcePn<B,
                                                      typename G::template Codim<0>::Entity,
                                                      typename G::ctype,
                                                      G::dimension,
                                                      DiscreteFunction<FvProductSpace<typename G::LeafGridView,
                                                                                      double,
                                                                                      rangeDim,
                                                                                      1>,
                                                                       typename Dune::XT::LA::
                                                                           Container<double,
                                                                                     XT::LA::default_sparse_backend>::
                                                                               VectorType>,
                                                      R,
                                                      rangeDim>>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d = G::dimension;
  static_assert(d == 3, "Only implemented for dimension 1.");
  typedef typename Hyperbolic::Problems::KineticEquation<
      typename Problems::
          PointSourcePn<B,
                        E,
                        D,
                        d,
                        DiscreteFunction<FvProductSpace<typename G::LeafGridView, double, rangeDim, 1>,
                                         typename Dune::XT::LA::Container<double,
                                                                          XT::LA::default_sparse_backend>::VectorType>,
                        R,
                        rangeDim>>
      ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::NonStationaryTestCase<G, ProblemType> BaseType;

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
