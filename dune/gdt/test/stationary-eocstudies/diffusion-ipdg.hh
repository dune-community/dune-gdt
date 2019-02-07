// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_DIFFUSION_IPDG_HH
#define DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_DIFFUSION_IPDG_HH

#include <dune/xt/la/container.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/filters/intersection.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/norms.hh>
#include <dune/gdt/operators/constant.hh>
#include <dune/gdt/operators/lincomb.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/oswald-interpolation.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \todo Add treatment of nonzero Dirichlet boundary values
 * \todo Add treatment of Neumann boundary values
 */
template <class G,
          LocalEllipticIpdgIntegrands::Method ipdg_method = LocalEllipticIpdgIntegrands::Method::swipdg,
          XT::LA::Backends la = XT::LA::Backends::istl_sparse>
class StationaryDiffusionIpdgEocStudy
  : public StationaryEocStudy<typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type,
                              1,
                              la>
{
protected:
  using BaseType =
      StationaryEocStudy<typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type, 1, la>;

  using BaseType::d;
  using typename BaseType::DF;
  using typename BaseType::E;
  using typename BaseType::GP;
  using typename BaseType::GV;
  using typename BaseType::I;
  using typename BaseType::M;
  using typename BaseType::O;
  using typename BaseType::S;
  using typename BaseType::V;

  using FF = XT::Functions::GridFunctionInterface<E>;
  using FT = XT::Functions::GridFunctionInterface<E, d, d>;

public:
  StationaryDiffusionIpdgEocStudy()
    : BaseType()
    , space_type_("")
  {}

protected:
  virtual const XT::Grid::BoundaryInfo<I>& boundary_info() const = 0;

  virtual const FF& diffusion_factor() const = 0;

  virtual const FT& diffusion_tensor() const = 0;

  virtual const FF& force() const = 0;

  virtual std::vector<std::string> norms() const override
  {
    auto nrms = BaseType::norms();
    nrms.push_back("eta_NC");
    return nrms;
  }

  virtual std::map<std::string, std::map<std::string, double>>
  compute(const size_t refinement_level,
          const std::vector<std::string>& actual_norms,
          const std::vector<std::pair<std::string, std::string>>& actual_estimates,
          const std::vector<std::string>& actual_quantities) override
  {
    auto& self = *this;
    // compute the quantities/norms/estimates we known about and remove them from the todos
    // store the data in current_data_, BaseType::compute will return that
    auto remaining_norms = actual_norms;
    auto remaining_estimates = actual_estimates;
    auto remaining_quantities = actual_quantities;
    for (auto norm_it = remaining_norms.begin(); norm_it != remaining_norms.end(); /*Do not increment here ...*/) {
      const auto norm_id = *norm_it;
      if (norm_id == "eta_NC") {
        norm_it = remaining_norms.erase(norm_it); // ... but rather here ...
        if (self.current_refinement_ != refinement_level)
          self.discretization_info(refinement_level);
        DUNE_THROW_IF(!self.current_space_, InvalidStateException, "");
        // compute current solution
        const auto& current_space = *self.current_space_;
        Timer timer;
        const auto solution = make_discrete_function(current_space, self.solve(current_space));
        // only set time if this did not happen in solve()
        if (self.current_data_["quantity"].count("time to solution (s)") == 0)
          self.current_data_["quantity"]["time to solution (s)"] = timer.elapsed();
        // compute estimate
        const auto h1_interpolation = apply_oswald_interpolation(solution, current_space, boundary_info());
        self.current_data_["norm"][norm_id] = elliptic_norm(
            current_space.grid_view(), diffusion_factor(), diffusion_tensor(), solution - h1_interpolation);
      } else
        ++norm_it; // ... or finally here.
    } // norms
    // let the Base compute the rest
    return BaseType::compute(refinement_level, remaining_norms, remaining_estimates, remaining_quantities);
  } // ... compute(...)

  virtual std::unique_ptr<S> make_space(const GP& current_grid) override
  {
    if (space_type_ == "fv")
      return std::make_unique<FiniteVolumeSpace<GV>>(current_grid.leaf_view());
    else if (space_type_.size() >= 4 && space_type_.substr(0, 4) == "dg_p") {
      const auto order = XT::Common::from_string<int>(space_type_.substr(4));
      return std::make_unique<DiscontinuousLagrangeSpace<GV>>(current_grid.leaf_view(), order);
    } else if (space_type_.size() >= 4 && space_type_.substr(0, 4) == "cg_p") {
      const auto order = XT::Common::from_string<int>(space_type_.substr(4));
      return std::make_unique<ContinuousLagrangeSpace<GV>>(current_grid.leaf_view(), order);
    } else {
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given, "space_type_ = " << space_type_);
      return nullptr;
    }
  } // ... make_space(...)

  virtual std::unique_ptr<O> make_residual_operator(const S& space) override
  {
    // define lhs operator (has to be a pointer to allow the residual operator to manage the memory in the end)
    auto lhs_op = std::make_unique<MatrixOperator<M, GV>>(make_matrix_operator<M>(
        space,
        (space_type_.size() >= 2 && space_type_.substr(0, 2) == "cg") ? Stencil::element
                                                                      : Stencil::element_and_intersection));
    lhs_op->append(
        LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(diffusion_factor(), diffusion_tensor())));
    lhs_op->append(LocalIntersectionIntegralBilinearForm<I>(LocalEllipticIpdgIntegrands::Inner<I, double, ipdg_method>(
                       diffusion_factor(), diffusion_tensor())),
                   {},
                   XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
    lhs_op->append(
        LocalIntersectionIntegralBilinearForm<I>(
            LocalEllipticIpdgIntegrands::DirichletBoundaryLhs<I, double, ipdg_method>(diffusion_factor(),
                                                                                      diffusion_tensor())),
        {},
        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(boundary_info(), new XT::Grid::DirichletBoundary()));
    // define rhs functional
    auto rhs_func = make_vector_functional<V>(space);
    rhs_func.append(LocalElementIntegralFunctional<E>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E>(), force())));
    // ... add Dirichlet here
    // (if we add something here, the oswald interpolation in compute() needs to be adapted accordingly!)
    // ... add Neumann here
    // assemble everything in one grid walk
    lhs_op->append(rhs_func);
    lhs_op->assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    // build residual operator
    auto residual_op = std::make_unique<ConstLincombOperator<M, GV>>(space, space);
    residual_op->add(lhs_op.release(), 1.);
    residual_op->add(new ConstantOperator<M, GV>(space, space, new V(std::move(rhs_func.vector()))), -1);
    return residual_op;
  } // ... make_residual_operator(...)

  std::string space_type_;
}; // class StationaryDiffusionIpdgEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_DIFFUSION_IPDG_HH
