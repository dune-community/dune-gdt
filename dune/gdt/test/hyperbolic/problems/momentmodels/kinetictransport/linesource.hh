// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_LINESOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_LINESOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/grid/view/periodic.hh>

#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

#include "../lebedevquadrature.hh"
#include "kinetictransportequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


// knots are the 20 linearly spaced knots between -1 and 1
// see https://de.wikipedia.org/wiki/Spline#B-Splines
template <size_t degree, size_t i, class D, class R>
struct BSpline
{
  static constexpr D h = 2. / 19.;

  static constexpr D t(const size_t j)
  {
    return -1. + j * h;
  }

  static constexpr R evaluate(const D& x)
  {
    return BSpline<degree - 1, i, D, R>::evaluate(x) * (x - t(i)) / (t(i + degree) - t(i))
           + BSpline<degree - 1, i + 1, D, R>::evaluate(x) * (t(i + degree + 1) - x) / (t(i + degree + 1) - t(i + 1));
  }
};

template <size_t i, class D, class R>
struct BSpline<0, i, D, R>
{
  static_assert(i < 19, "");

  static constexpr D h = 2. / 19.;

  static constexpr D t(const size_t j)
  {
    return -1. + j * h;
  }

  static constexpr R evaluate(const D& x)
  {
    return XT::Common::FloatCmp::ge(x, t(i)) && XT::Common::FloatCmp::le(x, t(i + 1)) ? 1 : 0;
  }
};

template <class BasisfunctionImp,
          class GridLayerImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class LineSourcePn : public KineticTransportEquation<BasisfunctionImp,
                                                     GridLayerImp,
                                                     EntityImp,
                                                     DomainFieldImp,
                                                     dimDomain,
                                                     U_,
                                                     RangeFieldImp,
                                                     dimRange,
                                                     dimDomain + 1>
{
  typedef KineticTransportEquation<BasisfunctionImp,
                                   GridLayerImp,
                                   EntityImp,
                                   DomainFieldImp,
                                   dimDomain,
                                   U_,
                                   RangeFieldImp,
                                   dimRange,
                                   dimDomain + 1>
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
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  LineSourcePn(const BasisfunctionType& basis_functions,
               const GridLayerType& grid_layer,
               const QuadratureType quadrature = QuadratureType(),
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, {1, 1}, grid_cfg, boundary_cfg, 1.)
  {
  }

  static std::string static_id()
  {
    return "modifiedlinesourcepn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-1 -1]";
    grid_config["upper_right"] = "[1 1]";
    grid_config["num_elements"] = "[5 5]";
    grid_config["overlap_size"] = "[1 1]";
    grid_config["num_quad_cells"] = "[10]";
    grid_config["quad_order"] = "50";
    return grid_config;
  }

  // sigma_a = 0, sigma_s = 0, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{0}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.5}),
                                  std::make_pair("t_end", std::vector<double>{0.05}),
                                  std::make_pair("num_segments", std::vector<double>{1., 1.})});
  }

  // Initial value of the kinetic equation is 1+p(||x||_2^2-0.1^2), where p is the spline of degree fourteen
  // on the twenty linearly-spaced knots over [-1,1] with coefficients [0, 0, 2.8309, 0, 0]
  // Thus the initial value for the moments is basis_integrated * (1+p(||x||_2^2-0.1^2))
  virtual InitialValueType* create_initial_values() const
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    RangeType basis_integrated = basis_functions_.integrated();
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = basis_integrated;
          ret *= 1 + 2.8309 * BSpline<14, 2, DomainFieldType, RangeFieldType>::evaluate(x.two_norm2() - 0.01);
          return ret;
        },
        15);
    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::grid_cfg_;
  using BaseType::basis_functions_;
  using BaseType::num_segments_;
  using BaseType::psi_vac_;
}; // class LineSourcePn<...>

template <class BasisfunctionType,
          class GridLayerType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class LineSourceMn : public LineSourcePn<BasisfunctionType,
                                         GridLayerType,
                                         EntityType,
                                         DomainFieldType,
                                         dimDomain,
                                         U_,
                                         RangeFieldType,
                                         dimRange>
{
  typedef LineSourcePn<BasisfunctionType,
                       GridLayerType,
                       EntityType,
                       DomainFieldType,
                       dimDomain,
                       U_,
                       RangeFieldType,
                       dimRange>
      BaseType;
  typedef LineSourceMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef GDT::EntropyBasedLocalFlux<BasisfunctionType,
                                     GridLayerType,
                                     EntityType,
                                     DomainFieldType,
                                     dimDomain,
                                     U_,
                                     RangeFieldType,
                                     dimRange,
                                     dimDomain + 1>
      ActualFluxType;
  using typename BaseType::QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  static QuadratureType default_quadrature(const XT::Common::Configuration& grid_cfg = default_grid_cfg())
  {
    size_t quad_order = grid_cfg.get("quad_order", 100);
    return LebedevQuadrature<DomainFieldType, true>::get(quad_order);
  }

  LineSourceMn(const BasisfunctionType& basis_functions,
               const GridLayerType& grid_layer,
               const QuadratureType& quadrature = default_quadrature(),
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "modifiedlinesourcemn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_, quadrature_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::quadrature_;
  using BaseType::grid_layer_;
}; // class LineSourceMn<...>


} // namespace KineticTransport
} // namespace Problems


template <class G,
          class R = double,
          size_t order = 6,
          class B = Hyperbolic::Problems::RealSphericalHarmonics<typename G::ctype, typename G::ctype, order, 2, true>>
class LineSourceTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               typename Hyperbolic::Problems::KineticEquation<
                                   typename Problems::KineticTransport::
                                       LineSourcePn<B,
                                                    typename G::LevelGridView,
                                                    typename G::template Codim<0>::Entity,
                                                    typename G::ctype,
                                                    G::dimension,
                                                    DiscreteFunction<FvProductSpace<typename G::LevelGridView,
                                                                                    double,
                                                                                    B::dimRange,
                                                                                    1>,
                                                                     typename Dune::XT::LA::
                                                                         Container<double,
                                                                                   XT::LA::default_sparse_backend>::
                                                                             VectorType>,
                                                    R,
                                                    B::dimRange>>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d = G::dimension;
  static_assert(d == 2, "Only implemented for dimension 2.");
  typedef typename Hyperbolic::Problems::KineticEquation<
      typename Problems::KineticTransport::
          LineSourcePn<B,
                       typename G::LevelGridView,
                       E,
                       D,
                       d,
                       DiscreteFunction<FvProductSpace<typename G::LevelGridView, double, B::dimRange, 1>,
                                        typename Dune::XT::LA::Container<double,
                                                                         XT::LA::default_sparse_backend>::VectorType>,
                       R,
                       B::dimRange>>
      ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;

  LineSourceTestCase(const size_t num_refs = 4, const double divide_t_end_by = 1.0)
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
        << "||  Testcase: LineSource Pn                                                                           ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1]^2                                                                                ||\n"
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
}; // class LineSourceTestCase

template <class G,
          class R = double,
          size_t order = 6,
          class B = Hyperbolic::Problems::RealSphericalHarmonics<typename G::ctype, typename G::ctype, order, 2, true>>
class LineSourceMnTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               typename Hyperbolic::Problems::KineticEquation<
                                   typename Problems::KineticTransport::LineSourceMn<
                                       typename XT::Grid::PeriodicGridView<typename G::LevelGridView>,
                                       B,
                                       typename G::template Codim<0>::Entity,
                                       typename G::ctype,
                                       G::dimension,
                                       DiscreteFunction<FvProductSpace<typename G::LevelGridView,
                                                                       double,
                                                                       B::dimRange,
                                                                       1>,
                                                        typename Dune::XT::LA::
                                                            Container<double,
                                                                      XT::LA::default_sparse_backend>::VectorType>,
                                       R,
                                       B::dimRange>>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d = G::dimension;
  static_assert(d == 2, "Only implemented for dimension 2.");
  typedef typename Hyperbolic::Problems::KineticEquation<typename Problems::KineticTransport::LineSourceMn<
      typename XT::Grid::PeriodicGridView<typename G::LevelGridView>,
      B,
      E,
      D,
      d,
      DiscreteFunction<FvProductSpace<typename G::LevelGridView, double, B::dimRange, 1>,
                       typename Dune::XT::LA::Container<double, XT::LA::default_sparse_backend>::VectorType>,
      R,
      B::dimRange>>
      ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;

  LineSourceMnTestCase(const size_t num_refs = 3, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
  {
  }

  virtual bool provides_exact_solution() const override final
  {
    return false;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+======================================================================================================+\n"
        << "|+====================================================================================================+|\n"
        << "||  Testcase: LineSource Mn                                                                           ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1]^2                                                                                ||\n"
        << "||  time = [0, 0.05]                                                                                  ||\n"
        << "||  flux = see http://dx.doi.org/10.1137/130934210 Section 6.5                                        ||\n"
        << "||  rhs = http://dx.doi.org/10.1137/130934210 Section 6.5                                             ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }
}; // class LineSourceTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_LINESOURCE_HH
