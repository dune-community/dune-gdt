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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_HH

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include "../kineticequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class BasisfunctionImp,
          class GridLayerImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          class RangeFieldImp,
          size_t rangeDim,
          size_t quadratureDim = domainDim>
class KineticTransportEquation : public KineticEquationImplementation<BasisfunctionImp,
                                                                      GridLayerImp,
                                                                      EntityImp,
                                                                      DomainFieldImp,
                                                                      domainDim,
                                                                      U_,
                                                                      RangeFieldImp,
                                                                      rangeDim>,
                                 public XT::Common::ParametricInterface
{
  typedef KineticTransportEquation ThisType;
  typedef KineticEquationImplementation<BasisfunctionImp,
                                        GridLayerImp,
                                        EntityImp,
                                        DomainFieldImp,
                                        domainDim,
                                        U_,
                                        RangeFieldImp,
                                        rangeDim>
      BaseType;

public:
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualFluxType;
  using typename BaseType::ActualRhsType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::RhsAffineFunctionType;

  typedef Dune::QuadratureRule<DomainFieldType, quadratureDim> QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  static QuadratureType default_quadrature(const XT::Common::Configuration& grid_cfg = default_grid_cfg())
  {
    std::vector<int> num_quad_cells = grid_cfg.get("num_quad_cells", std::vector<int>{2, 2, 2});
    int quad_order = grid_cfg.get("quad_order", 20);
    // quadrature that consists of a Gauss-Legendre quadrature on each cell of the velocity grid
    QuadratureType quadrature;
    Dune::FieldVector<double, quadratureDim> lower_left(-1);
    Dune::FieldVector<double, quadratureDim> upper_right(1);
    std::array<int, quadratureDim> s;
    for (size_t ii = 0; ii < quadratureDim; ++ii)
      s[ii] = num_quad_cells[ii];
    typedef typename Dune::YaspGrid<quadratureDim, Dune::EquidistantOffsetCoordinates<double, quadratureDim>> GridType;
    GridType velocity_grid(lower_left, upper_right, s);
    const auto velocity_grid_view = velocity_grid.leafGridView();
    for (const auto& entity : elements(velocity_grid_view)) {
      const auto local_quadrature = Dune::QuadratureRules<DomainFieldType, quadratureDim>::rule(
          entity.type(), quad_order, Dune::QuadratureType::GaussLegendre);
      for (const auto& quad_point : local_quadrature) {
        quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, quadratureDim>(
            entity.geometry().global(quad_point.position()),
            quad_point.weight() * entity.geometry().integrationElement(quad_point.position())));
      }
    }
    return quadrature;
  }

  KineticTransportEquation(const BasisfunctionType& basis_functions,
                           const GridLayerType& grid_layer,
                           const QuadratureType quadrature = QuadratureType(),
                           const FieldVector<size_t, dimDomain> num_segments = {1, 1, 1},
                           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
                           const RangeFieldType psi_vac = 5e-9,
                           const XT::Common::ParameterType& parameter_type = XT::Common::ParameterType())
    : BaseType(basis_functions, grid_layer)
    , num_segments_(num_segments)
    , grid_cfg_(grid_cfg)
    , boundary_cfg_(boundary_cfg)
    , psi_vac_(psi_vac)
    , quadrature_(quadrature)
    , parameter_type_(parameter_type)
  {
    if (quadrature_.empty())
      quadrature_ = default_quadrature(grid_cfg);
    if (parameter_type_.empty())
      parameter_type_ = XT::Common::ParameterType({std::make_pair("sigma_a", get_num_regions(num_segments)),
                                                   std::make_pair("sigma_s", get_num_regions(num_segments)),
                                                   std::make_pair("Q", get_num_regions(num_segments)),
                                                   std::make_pair("CFL", 1),
                                                   std::make_pair("t_end", 1)});
  }

  virtual ~KineticTransportEquation()
  {
  }

  virtual bool is_parametric() const override
  {
    return true;
  }

  virtual const XT::Common::ParameterType& parameter_type() const override
  {
    return parameter_type_;
  }

  virtual XT::Common::Parameter parameters() const = 0;

  using XT::Common::ParametricInterface::parse_parameter;

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  virtual FluxType* create_flux() const override
  {
    auto A = basis_functions_.mass_matrix_with_v();
    //    std::cout << XT::Common::to_string(A) << std::endl;
    auto M_inv = basis_functions_.mass_matrix_inverse();
    for (size_t dd = 0; dd < dimDomain; ++dd)
      A[dd].rightmultiply(M_inv);
    return new ActualFluxType(A, typename ActualFluxType::RangeType(0));
  }

  // RHS is (sigma_s/vol*G - sigma_t * I)u + Q<b>,
  // where sigma_t = sigma_s + sigma_a, G = <b><b>^T M^{-1} = <b>*c^T and
  // vol = <1> is the volume of the integration domain.
  virtual RhsType* create_rhs() const override
  {
    const auto param = parse_parameter(parameters());
    const size_t num_regions = get_num_regions(num_segments_);
    auto sigma_a = param.get("sigma_a");
    auto sigma_s = param.get("sigma_s");
    auto Q = param.get("Q");
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    auto sigma_t = sigma_a;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_t[ii] += sigma_s[ii];
    const RangeType basis_integrated = basis_functions_.integrated();
    const auto M_inv = basis_functions_.mass_matrix_inverse();
    RangeType c(0);
    M_inv.mtv(basis_integrated, c);
    MatrixType I(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType G(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        G[rr][cc] = basis_integrated[rr] * c[cc];
    const auto vol = unit_ball_volume();

    std::vector<RhsAffineFunctionType> affine_functions;
    for (size_t ii = 0; ii < num_regions; ++ii) {
      MatrixType G_scaled = G;
      G_scaled *= sigma_s[ii] / vol;
      MatrixType I_scaled = I;
      I_scaled *= sigma_t[ii];
      MatrixType A = G_scaled;
      A -= I_scaled;
      RangeType b = basis_integrated;
      b *= Q[ii];
      affine_functions.emplace_back(A, b, "rhs");
    } // ii
    return new ActualRhsType(lower_left, upper_right, num_segments_, affine_functions);
  } // ... create_rhs(...)

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is basis_integrated * psi_vac.
  virtual InitialValueType* create_initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    RangeType value = basis_functions_.integrated();
    value *= psi_vac_;
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    const size_t num_regions = get_num_regions(num_segments_);
    for (size_t ii = 0; ii < num_regions; ++ii)
      initial_vals.emplace_back([=](const DomainType&, const XT::Common::Parameter&) { return value; }, 0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

  // Use a constant vacuum concentration basis_integrated * psi_vac as boundary value
  virtual BoundaryValueType* create_boundary_values() const override
  {
    RangeType value = basis_functions_.integrated();
    value *= psi_vac_;
    return new ActualBoundaryValueType([=](const DomainType&, const XT::Common::Parameter&) { return value; }, 0);
  } // ... create_boundary_values()

  virtual RangeFieldType CFL() const override
  {
    return parameters().get("CFL")[0];
  }

  virtual RangeFieldType t_end() const override
  {
    return parameters().get("t_end")[0];
  }

  virtual XT::Common::Configuration grid_config() const override
  {
    return grid_cfg_;
  }

  virtual XT::Common::Configuration boundary_config() const override
  {
    return boundary_cfg_;
  }

  static std::string static_id()
  {
    return "kinetictransportequation";
  }

  virtual const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

protected:
  static size_t get_num_regions(const FieldVector<size_t, dimDomain>& num_segments)
  {
    return std::accumulate(num_segments.begin(), num_segments.end(), 1, [](auto a, auto b) { return a * b; });
  }

  static RangeFieldType unit_ball_volume()
  {
    if (dimDomain == 1)
      return 2;
    else if (dimDomain == 2)
      return 2 * M_PI;
    else if (dimDomain == 3)
      return 4 * M_PI;
    else {
      DUNE_THROW(NotImplemented, "");
      return 0;
    }
  }

  using BaseType::basis_functions_;
  const FieldVector<size_t, dimDomain> num_segments_;
  const XT::Common::Configuration grid_cfg_;
  const XT::Common::Configuration boundary_cfg_;
  const RangeFieldType psi_vac_;
  QuadratureType quadrature_;
  XT::Common::ParameterType parameter_type_;
}; // class KineticTransportEquation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_HH
