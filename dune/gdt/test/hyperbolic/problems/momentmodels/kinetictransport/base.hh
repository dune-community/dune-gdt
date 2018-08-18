// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH

#include <dune/grid/common/partitionset.hh>

#include "../kineticequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class BasisfunctionImp, class GridLayerImp, class U_, size_t quadratureDim = BasisfunctionImp::dimDomain>
class KineticTransportEquation : public KineticEquationImplementationInterface<BasisfunctionImp, GridLayerImp, U_>,
                                 public XT::Common::ParametricInterface
{
  typedef KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_, quadratureDim> ThisType;
  typedef KineticEquationImplementationInterface<BasisfunctionImp, GridLayerImp, U_> BaseType;

public:
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
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
    Dune::FieldVector<double, quadratureDim> lower_left(-1);
    Dune::FieldVector<double, quadratureDim> upper_right(1);
    std::array<int, quadratureDim> s;
    for (size_t ii = 0; ii < quadratureDim; ++ii)
      s[ii] = num_quad_cells[ii];
    typedef typename Dune::YaspGrid<quadratureDim, Dune::EquidistantOffsetCoordinates<double, quadratureDim>> GridType;
    GridType velocity_grid(lower_left, upper_right, s);
    const auto velocity_grid_view = velocity_grid.leafGridView();
    // In an MPI parallel environment, each rank will only have parts of the velocity grid,
    // but we need the full quadrature on each rank, so we need to communicate the missing parts
    std::vector<double> coords;
    std::vector<double> weights;
    // first, collect the quadrature points on each rank
    for (const auto& entity : elements(velocity_grid_view, Dune::Partitions::Interior())) {
      const auto local_quadrature = Dune::QuadratureRules<DomainFieldType, quadratureDim>::rule(
          entity.type(), quad_order, Dune::QuadratureType::GaussLobatto);
      for (const auto& quad_point : local_quadrature) {
        const auto pos = entity.geometry().global(quad_point.position());
        for (size_t ii = 0; ii < quadratureDim; ++ii)
          coords.push_back(pos[ii]);
        weights.push_back(quad_point.weight() * entity.geometry().integrationElement(quad_point.position()));
      }
    }
    // communicate the quadrature points to all ranks
    const auto& comm = velocity_grid.comm();
    assert(coords.size() < std::numeric_limits<int>::max() && weights.size() < std::numeric_limits<int>::max());
    int coords_size = static_cast<int>(coords.size());
    int weights_size = static_cast<int>(weights.size());
    assert(comm.size() >= 0 && coords_size >= 0 && weights_size >= 0);
    size_t comm_size = static_cast<size_t>(comm.size());
    std::vector<int> coords_sizes(comm_size);
    std::vector<int> weights_sizes(comm_size);
    std::vector<int> coords_displacements(comm_size, 0);
    std::vector<int> weights_displacements(comm_size, 0);
    comm.allgather(&coords_size, 1, coords_sizes.data());
    comm.allgather(&weights_size, 1, weights_sizes.data());
    for (size_t ii = 1; ii < comm_size; ++ii) {
      coords_displacements[ii] = coords_displacements[ii - 1] + coords_sizes[ii - 1];
      weights_displacements[ii] = weights_displacements[ii - 1] + weights_sizes[ii - 1];
    }
    size_t all_coords_size = static_cast<size_t>(comm.sum(coords_size));
    size_t all_weights_size = static_cast<size_t>(comm.sum(weights_size));
    std::vector<double> all_coords(all_coords_size);
    std::vector<double> all_weights(all_weights_size);
    comm.allgatherv(coords.data(), coords_size, all_coords.data(), coords_sizes.data(), coords_displacements.data());
    comm.allgatherv(
        weights.data(), weights_size, all_weights.data(), weights_sizes.data(), weights_displacements.data());
    // now assemble the quadrature
    QuadratureType quadrature;
    for (size_t ii = 0; ii < all_weights.size(); ++ii) {
      XT::Common::FieldVector<DomainFieldType, quadratureDim> position;
      for (size_t jj = 0; jj < quadratureDim; ++jj)
        position[jj] = all_coords[ii * quadratureDim + jj];
      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, quadratureDim>(position, all_weights[ii]));
    }
    return quadrature;
  }

  KineticTransportEquation(const BasisfunctionType& basis_functions,
                           const GridLayerType& grid_layer,
                           const QuadratureType quadrature = QuadratureType(),
                           DynamicVector<size_t> num_segments = {1, 1, 1},
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
    // by default, initially num_segments has size 3. If dimDomain < 3 we thus need to resize
    num_segments_.resize(dimDomain);
    if (quadrature_.empty())
      quadrature_ = default_quadrature(grid_cfg);
    if (parameter_type_.empty())
      parameter_type_ = XT::Common::ParameterType({std::make_pair("sigma_a", get_num_regions(num_segments_)),
                                                   std::make_pair("sigma_s", get_num_regions(num_segments_)),
                                                   std::make_pair("Q", get_num_regions(num_segments_)),
                                                   std::make_pair("CFL", 1),
                                                   std::make_pair("t_end", 1)});
  }

  virtual ~KineticTransportEquation() override
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
    // calculate B row-wise by solving M^{T} A^T = B^T column-wise
    auto A = basis_functions_.mass_matrix_with_v();
    const auto M_T = basis_functions_.mass_matrix(); // mass matrix is symmetric
    // solve
    DynamicVector<RangeFieldType> tmp_row(M_T.N(), 0.);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t ii = 0; ii < M_T.N(); ++ii) {
        M_T.solve(tmp_row, A[dd][ii]);
        A[dd][ii] = tmp_row;
      }
    }
    return new ActualFluxType(A, typename ActualFluxType::RangeType(0));
  }

  // RHS is (sigma_s/vol*G - sigma_t * I)u + Q<b>,
  // where sigma_t = sigma_s + sigma_a, G = <b><b>^T M^{-1} = <b>*c^T and
  // vol = <1> is the volume of the integration domain.
  virtual RhsType* create_rhs() const override
  {
    const auto param = parse_parameter(parameters());
    const size_t num_regions = get_num_regions(num_segments_);
    const auto sigma_a = param.get("sigma_a");
    const auto sigma_s = param.get("sigma_s");
    const auto Q = param.get("Q");
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    auto sigma_t = sigma_a;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_t[ii] += sigma_s[ii];
    const RangeType basis_integrated = basis_functions_.integrated();
    // calculate c = M^{-T} <b>
    const auto M_T = basis_functions_.mass_matrix(); // mass matrix is symmetric
    RangeType c(0.);
    M_T.solve(c, basis_integrated);
    MatrixType I(dimRange, dimRange, 0.);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType G(dimRange, dimRange, 0.);
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
      affine_functions.emplace_back(A, b, true, "rhs");
    } // ii
    return new ActualRhsType(lower_left, upper_right, num_segments_, affine_functions);
  } // ... create_rhs(...)

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is basis_integrated * psi_vac.
  virtual InitialValueType* create_initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    RangeType value = basis_functions_.integrated() * psi_vac_;
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    const size_t num_regions = get_num_regions(num_segments_);
    for (size_t ii = 0; ii < num_regions; ++ii)
      initial_vals.emplace_back([=](const DomainType&, const XT::Common::Parameter&) { return value; }, 0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

  // Use a constant vacuum concentration basis_integrated * psi_vac as boundary value
  virtual BoundaryValueType* create_boundary_values() const override
  {
    RangeType value = basis_functions_.integrated() * psi_vac_;
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

  virtual const RangeFieldType psi_vac() const
  {
    return psi_vac_;
  }

protected:
  static size_t get_num_regions(const DynamicVector<size_t>& num_segments)
  {
    return std::accumulate(num_segments.begin(), num_segments.end(), size_t(1), [](auto a, auto b) { return a * b; });
  }

  static RangeFieldType unit_ball_volume()
  {
    if (BasisfunctionType::dimDomain == 1)
      return 2;
    else if (BasisfunctionType::dimDomain == 2)
      return 2 * M_PI;
    else if (BasisfunctionType::dimDomain == 3)
      return 4 * M_PI;
    else {
      DUNE_THROW(NotImplemented, "");
      return 0;
    }
  }

  using BaseType::basis_functions_;
  DynamicVector<size_t> num_segments_;
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

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH
