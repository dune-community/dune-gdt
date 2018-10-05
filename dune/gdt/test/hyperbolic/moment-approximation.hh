// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_MOMENT_APPROXIMATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_MOMENT_APPROXIMATION_HH

#include <chrono>
#include <fstream>

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/operators/fv/quadrature.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/kineticequation.hh>

#include "pn-discretization.hh"

template <class TestCaseType>
struct MomentApproximation
{
  static void run(size_t num_quad_refinements = TestCaseType::RealizabilityLimiterChooserType::num_quad_refinements,
                  size_t quad_order = TestCaseType::RealizabilityLimiterChooserType::quad_order,
                  std::string filename = "")
  {
    using namespace Dune;
    using namespace Dune::GDT;

    //******************* get typedefs and constants from ProblemType **********************//
    using BasisfunctionType = typename TestCaseType::BasisfunctionType;
    using DiscreteFunctionType = typename TestCaseType::DiscreteFunctionType;
    using GridType = typename TestCaseType::GridType;
    using SpaceType = typename TestCaseType::SpaceType;
    using GridLayerType = typename TestCaseType::GridLayerType;
    using ProblemType = typename TestCaseType::ProblemType;
    using EquationType = Hyperbolic::Problems::KineticEquation<ProblemType>;
    using DomainType = typename EquationType::DomainType;
    using RangeType = typename EquationType::RangeType;
    using RangeFieldType = typename EquationType::RangeFieldType;
    static constexpr size_t dimDomain = BasisfunctionType::dimDomain;
    static constexpr size_t dimRange = BasisfunctionType::dimRange;

    //******************* create grid and FV space ***************************************
    auto grid_config = EquationType::default_grid_cfg();
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GridLayerType grid_layer(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_layer);

    //******************* create EquationType object ***************************************
    std::shared_ptr<const BasisfunctionType> basis_functions =
        std::make_shared<const BasisfunctionType>(quad_order, num_quad_refinements);
    const std::unique_ptr<ProblemType> problem_imp =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_layer, grid_config);
    const EquationType problem(*problem_imp);

    // ***************** project initial values to discrete function *********************
    // Heaviside
    // const auto psi = [](const DomainType& v) { return v[0] < 0 ? RangeFieldType(5e-9) : RangeFieldType(1); };
    // Gauss
    // const RangeFieldType sigma = 0.5;
    // const RangeFieldType mu = 0.;
    // const auto psi = [sigma, mu](const DomainType& v){ return 1./std::sqrt(2. * M_PI * std::pow(sigma, 2)) *
    // std::exp(std::pow(v[0]-mu, 2)/(-2 *std::pow(sigma, 2))); };
    // Gauss on sphere
    // const RangeFieldType sigma = 0.5;
    // const RangeFieldType mu = 1.;
    // const auto psi = [sigma, mu](const DomainType& v){ return 1./std::sqrt(2. * M_PI * std::pow(sigma, 2)) *
    // std::exp(std::pow(v[0]-mu, 2)/(-2 *std::pow(sigma, 2))); };
    // Square on sphere
    const auto psi = [](const DomainType& v) {
      if (v[0] > 0 && std::abs(v[1]) < 0.5 && std::abs(v[2]) < 0.5)
        return RangeFieldType(1);
      else
        return RangeFieldType(5e-9);
    };
    const RangeType u = basis_functions->get_moment_vector(psi);
    // const RangeType u{2.82081e-09,  1.47519e-13,  1.47552e-13,  -1.86436e-16, 1.32958e-19,  -1.8985e-13,
    // -4.30881e-14,
    //                  2.13323e-16,  7.47148e-14,  -2.13786e-14, -2.80844e-19, 1.18799e-13,  -5.60185e-14,
    //                  -8.48068e-17,
    //                  -1.06901e-13, -6.78253e-17, 5.92712e-20,  3.28499e-14,  4.21767e-19,  -1.24544e-14, 7.72642e-14,
    //                  -6.15351e-17, 8.53369e-14,  9.48222e-17,  6.43495e-16,  -2.62515e-15, 9.43676e-20, -2.99768e-14,
    //                  -3.62141e-19, -4.25964e-14, -3.82627e-14, 1.00791e-16,  -3.56862e-14, -6.94389e-17,
    //                  -1.09698e-15,
    //                  -6.16735e-18, -4.9294e-19,  4.4688e-15,   -5.17249e-19, 1.7311e-14,   8.12798e-20,  3.34336e-14,
    //                  -5.09742e-16, -4.25455e-17, 3.16988e-16,  1.80488e-17,  1.17459e-15,  9.62318e-18,  2.38981e-16,
    //                  -8.14399e-16, 9.28639e-19,  -4.8444e-15,  1.01707e-18,  -6.29709e-15, 2.89504e-19, -5.97879e-15,
    //                  7.88486e-15,  -1.78209e-17, 5.66809e-15,  1.45894e-17,  -1.08493e-15, -8.51736e-18,
    //                  -4.39289e-16,
    //                  -2.53324e-19, 3.32864e-19,  1.47777e-15,  -9.93208e-19, 4.01022e-15,  -1.27832e-18, 2.71052e-15,
    //                  -4.61338e-19, -2.48712e-15, 6.90215e-16,  2.37436e-17,  7.92837e-16,  -1.41213e-17, 3.56008e-16,
    //                  5.12667e-18,  2.46417e-16,  5.24185e-19,  -1.44969e-17};

    using MnFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, DiscreteFunctionType>;
    MnFluxType mn_flux(*basis_functions, grid_layer);
    auto local_mn_flux = mn_flux.derived_local_function(*(grid_layer.template begin<0>()));
    const auto alpha = local_mn_flux->get_alpha(DomainType(0), u, {"boundary", {0}}, false).first;

    std::string filename_mn = filename;
    filename += "_" + basis_functions->short_id();
    std::ofstream mn_file(filename + "_m" + Dune::XT::Common::to_string(dimRange) + ".txt");
    std::ofstream pn_file(filename + "_p" + Dune::XT::Common::to_string(dimRange) + ".txt");
    std::ofstream mn_errors_file(filename + "_mn_errors.txt", std::ios_base::app);
    std::ofstream pn_errors_file(filename + "_pn_errors.txt", std::ios_base::app);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      mn_file << "v" << dd << " ";
      pn_file << "v" << dd << " ";
    }
    mn_file << "psi" << std::endl;
    pn_file << "psi" << std::endl;
    RangeFieldType l1error_mn(0), l2error_mn(0), linferror_mn(0);
    RangeFieldType l1error_pn(0), l2error_pn(0), linferror_pn(0);
    const auto mass_matrix = basis_functions->mass_matrix();
    RangeType pn_coeffs;
    mass_matrix.solve(pn_coeffs, u);
    const auto quadrature = basis_functions->quadratures().merged();
    for (auto it = quadrature.begin(); it != quadrature.end(); ++it) {
      const auto& quad_point = *it;
      const auto v = quad_point.position();
      const auto weight = quad_point.weight();
      const auto basis = basis_functions->evaluate(v, it.first_index());
      const auto psi_mn = std::exp(alpha * basis);
      const auto psi_pn = pn_coeffs * basis;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        mn_file << XT::Common::to_string(v[dd], 15) << " ";
        pn_file << XT::Common::to_string(v[dd], 15) << " ";
      }
      mn_file << XT::Common::to_string(psi_mn, 15) << std::endl;
      pn_file << XT::Common::to_string(psi_pn, 15) << std::endl;
      l1error_mn += std::abs(psi_mn - psi(v)) * weight;
      l2error_mn += std::pow(psi_mn - psi(v), 2) * weight;
      linferror_mn = std::max(std::abs(psi_mn - psi(v)), linferror_mn);
      l1error_pn += std::abs(psi_pn - psi(v)) * weight;
      l2error_pn += std::pow(psi_pn - psi(v), 2) * weight;
      linferror_pn = std::max(std::abs(psi_pn - psi(v)), linferror_pn);
    }
    mn_errors_file << dimRange << " " << XT::Common::to_string(l1error_mn, 15) << " "
                   << XT::Common::to_string(l2error_mn, 15) << " " << XT::Common::to_string(linferror_mn, 15)
                   << std::endl;
    pn_errors_file << dimRange << " " << XT::Common::to_string(l1error_pn, 15) << " "
                   << XT::Common::to_string(l2error_pn, 15) << " " << XT::Common::to_string(linferror_pn, 15)
                   << std::endl;
  }
};


#endif // DUNE_GDT_TEST_HYPERBOLIC_MOMENT_APPROXIMATION_HH
