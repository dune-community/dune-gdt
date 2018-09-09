// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH

#include <dune/xt/common/string.hh>

#include <dune/xt/grid/information.hh>

#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/operators/fv/reconstruction/linear.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class PlaneSourcePn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>
{
  using BaseType = KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualDirichletBoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::FluxType;
  using BaseType::default_boundary_cfg;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  PlaneSourcePn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, {1}, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "planesourcepn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-1.2]";
    grid_config["upper_right"] = "[1.2]";
    grid_config["num_elements"] = "[240]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  // sigma_a = 0, sigma_s = 1, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{1}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.49}),
                                  std::make_pair("t_end", std::vector<double>{1.0})});
  }

  // Initial value of the kinetic equation is psi_vac + delta(x).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + delta(x))
  virtual InitialValueType* create_initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    assert(grid_layer_.size(0) >= 0 && grid_layer_.overlapSize(0) >= 0
           && grid_layer_.overlapSize(0) < grid_layer_.size(0));
    size_t num_elements = static_cast<size_t>(grid_layer_.size(0));
    if (grid_layer_.comm().size() > 1) {
      num_elements -= static_cast<size_t>(grid_layer_.overlapSize(0));
      num_elements = grid_layer_.comm().sum(num_elements);
    }
    if (num_elements % 2)
      DUNE_THROW(Dune::NotImplemented, "An even number of grid cells is needed for this test!");
    const RangeFieldType len_domain = upper_right[0] - lower_left[0];
    const RangeFieldType vol_entity = len_domain / num_elements;
    RangeType basis_integrated = basis_functions_.integrated();
    const RangeFieldType domain_center = lower_left[0] + len_domain / 2;
    // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = basis_integrated;
          if (XT::Common::FloatCmp::ge(x[0], domain_center - vol_entity)
              && XT::Common::FloatCmp::le(x[0], domain_center + vol_entity))
            ret *= psi_vac_ + 1. / (2. * vol_entity);
          else
            ret *= psi_vac_;
          return ret;
        },
        0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

  // return exact solution if there is no rhs (i.e. sigma_s = sigma_a = Q = 0) and psi_vac = 0
  std::unique_ptr<InitialValueType> exact_solution_without_rhs() const
  {
    const std::unique_ptr<FluxType> flux(BaseType::create_flux());
    const auto jacobian =
        flux->local_function(grid_layer_.template begin<0>())->partial_u_col(0, DomainType(0.), RangeType(0.));
    using MatrixType = typename FluxType::PartialURangeType;
    using EigenSolverType = Dune::XT::LA::EigenSolver<MatrixType>;
    static auto eigensolver_options = Dune::GDT::internal::hyperbolic_default_eigensolver_options<MatrixType>();
    const auto eigensolver = EigenSolverType(jacobian, &eigensolver_options);
    const auto eigenvectors = eigensolver.real_eigenvectors();
    const auto eigenvalues = eigensolver.real_eigenvalues();
    const auto eigenvectors_inv = eigensolver.real_eigenvectors_inverse();
    XT::Grid::Dimensions<GridLayerType> dimensions(grid_layer_);
    RangeFieldType dx = dimensions.entity_width.max();
    RangeType initial_val = basis_functions_.integrated() * 1. / (2. * dx);
    RangeType char_initial_val;
    eigenvectors_inv.mv(initial_val, char_initial_val);
    auto lambda = [=](const FieldVector<double, 1>& x, const XT::Common::Parameter& param) {
      const auto t = param.get("t")[0];
      RangeType ret(0.);
      for (size_t ii = 0; ii < dimRange; ++ii) {
        auto x_initial_ii = x[0] - eigenvalues[ii] * t;
        if (Dune::XT::Common::FloatCmp::le(x_initial_ii, dx, 1e-8, 1e-8)
            && Dune::XT::Common::FloatCmp::ge(x_initial_ii, -dx, 1e-8, 1e-8))
          ret[ii] = char_initial_val[ii];
        else
          ret[ii] = 0.;
      }
      const auto ret_copy = ret;
      eigenvectors.mv(ret_copy, ret);
      return ret;
    };
    return std::make_unique<XT::Functions::GlobalLambdaFunction<typename GridLayerType::template Codim<0>::Entity,
                                                                DomainFieldType,
                                                                dimDomain,
                                                                RangeFieldType,
                                                                dimRange,
                                                                1>>(lambda, 10);
  }

protected:
  using BaseType::grid_cfg_;
  using BaseType::grid_layer_;
  using BaseType::basis_functions_;
  using BaseType::num_segments_;
  using BaseType::psi_vac_;
}; // class PlaneSourcePn<...>


template <class BasisfunctionType, class GridLayerImp, class U_>
class PlaneSourceMn : public PlaneSourcePn<BasisfunctionType, GridLayerImp, U_>
{
  using BaseType = PlaneSourcePn<BasisfunctionType, GridLayerImp, U_>;
  using ThisType = PlaneSourceMn;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  using ActualFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerImp, U_>;
  using typename BaseType::GridLayerType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  PlaneSourceMn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "planesourcemn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
}; // class PlaneSourceMn<...>


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
