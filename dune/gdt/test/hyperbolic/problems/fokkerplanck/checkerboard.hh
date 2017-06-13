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
class CheckerboardPn : public KineticTransportEquation<BasisfunctionImp,
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

  CheckerboardPn(const BasisfunctionType& basis_functions,
                 const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                 const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "checkerboardpn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0 0 0]";
    grid_config["upper_right"] = "[7 7 7]";
    grid_config["num_elements"] = "[7 7 7]";
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
                                  std::make_pair("t_end", std::vector<double>{3.2}),
                                  std::make_pair("num_segments", std::vector<double>{1., 1., 1.})});
  }

  // sigma_a = 0, sigma_s = 1, Q = 0 in scattering regions
  // sigma_a = 10, sigma_s = 0, Q = 0 in absorbing regions
  // sigma_a = 0, sigma_s = 1, Q = 1 for center cube
  static ConfigType create_rhs_config(const ConfigType grid_config,
                                      const Dune::QuadratureRule<double, 3>& quadrature,
                                      const SphericalTriangulation<double>& poly)
  {
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[7 7 7]";
    rhs_config["name"] = DefaultRHSType::static_id();

    RangeType Q;
    RangeFieldType sigma_s, sigma_t;

    for (size_t plane = 0; plane < 7; ++plane) {
      for (size_t row = 0; row < 7; ++row) {
        for (size_t col = 0; col < 7; ++col) {
          if (plane == 3 && row == 3 && col == 3) { // center
            Q = basis_integrated;
            sigma_s = 1;
            sigma_t = 1;
          } else if (is_absorbing(plane, row, col)) { // absorbing regions
            Q *= 0;
            sigma_s = 0;
            sigma_t = 10;
          } else { // scattering regions (without center)
            Q *= 0;
            sigma_s = 1;
            sigma_t = 1;
          }

          MatrixType A(1);
          A *= sigma_s / (4 * M_PI);
          for (size_t nn = 0; nn < dimRange; ++nn) {
            A[nn] *= basis_integrated[nn];
            A[nn][nn] -= sigma_t;
          }
          size_t number = 49 * plane + 7 * row + col;
          rhs_config["A." + Dune::XT::Common::to_string(number)] = Dune::XT::Common::to_string(A, precision);
          rhs_config["b." + Dune::XT::Common::to_string(number)] = Dune::XT::Common::to_string(Q, precision);
        } // col
      } // row
    } // plane

    return rhs_config;
  } // ... create_rhs_config()

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
        << "||  Testcase: PointSource Pn                                                                ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [-0.5, 0.5]^3                                                                            ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end()) + "]                                ||\n"
        << "||  flux = see http://dx.doi.org/10.1137/130934210 Section 6.5                                        ||\n"
        << "||  rhs = http://dx.doi.org/10.1137/130934210 Section 6.5                                             ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const ProblemType problem_;
}; // class SourceBeamTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
