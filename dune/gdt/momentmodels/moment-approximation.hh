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

#include <dune/xt/common/string.hh>

#include <dune/xt/data/quadratures/fekete.hh>

#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/momentmodels/entropyflux.hh>

namespace Dune {
namespace GDT {


bool is_empty(const std::string filename)
{
  std::ifstream file(filename);
  return !file.good() || file.peek() == std::ifstream::traits_type::eof();
}

template <class MomentBasisType, class DiscreteFunctionType>
struct MomentApproximation
{
  using SpaceType = typename DiscreteFunctionType::SpaceType;
  using GridViewType = typename SpaceType::GridViewType;
  using GridType = typename GridViewType::Grid;
  using DomainType = typename MomentBasisType::DomainType;
  using RangeFieldType = typename DiscreteFunctionType::RangeFieldType;
  using DynamicRangeType = typename MomentBasisType::DynamicRangeType;
  using QuadraturesType = typename MomentBasisType::QuadraturesType;
  static const size_t dimDomain = MomentBasisType::dimDomain;
  static const size_t dimRange = MomentBasisType::dimRange;

  // Chooses the quadrature that is used to calculate l1, l2, linf error and do the visualization
  // In 1d, we use a Gauss-Lobatto-Quadrature of order 197 (100 quadrature points) on each interval. For the
  // fullmoments, we fix the number of intervals to 50. Thus we have a different quadrature for each model, but the high
  // order of the quadratures should ensure that the quadrature error is negligible compared to the model error.
  // Choosing a fixed quadrature for all models would not allow to align the quadrature to the interval boundaries,
  // which potentially causes much higher errors especially in linf norm. In 3d, all spherical triangles are refinements
  // of the initial octants, so we can use the same quadrature (on the finest triangulation) for all models. Note that
  // this quadrature is not used for the optimization problems in the entropy-based models. For the optimization
  // problems, we use the following model-adapted quadratures: PMM_n and HFM_n (1d): A Gauss-Lobatto quadrature of order
  // 15 on each interval. M_N (1d): A Gauss-Lobatto quadrature of order 197 on each of the intervals [-1,0] and [0,1].
  // PMM_n and HFM_n (3d): A Fekete quadrature of order 9 on each spherical triangle (order 15 for PMM_32 and HFM_6).
  // M_N (3d): A product quadrature of order 2N+8 on the octants of the sphere.
  template <size_t domainDim = dimDomain, bool anything = true>
  struct QuadratureHelper;

  template <bool anything>
  struct QuadratureHelper<3, anything>
  {
    static QuadraturesType get(const int quadrature_refinements, const int /*additional_refs*/)
    {
      QuadratureRule<RangeFieldType, 2> reference_quadrature_rule = XT::Data::FeketeQuadrature<RangeFieldType>::get(1);
      int quad_refs = quadrature_refinements;
      const size_t num_refs = MomentBasisType::num_refinements;
      const int additional_refs = quad_refs - static_cast<int>(num_refs);
      if (additional_refs < 0)
        DUNE_THROW(InvalidStateException,
                   "The number of refinements for the quadrature has to be greater or equal than the number of "
                   "refinements the basisfunctions use!");
      SphericalTriangulation<RangeFieldType> triangulation(num_refs);
      return triangulation.quadrature_rules(static_cast<size_t>(additional_refs), reference_quadrature_rule);
    }
  };

  template <bool anything>
  struct QuadratureHelper<1, anything>
  {
    static QuadraturesType get(const int num_quadrature_intervals, const int additional_refs)
    {
      return MomentBasisType::gauss_lobatto_quadratures(num_quadrature_intervals, 197, additional_refs);
    }
  };

  // use block-wise solver for 3d partial moments, sparse solver for full matrix otherwise
  template <bool partialmoments = is_3d_partial_moment_basis<MomentBasisType>::value, bool anything = true>
  struct SolverHelper;

  template <bool anything>
  struct SolverHelper<true, anything>
  {
    static void solve(const MomentBasisType& basis_functions, const DynamicRangeType& u, DynamicRangeType& pn_coeffs)
    {
      const auto mass_matrix = basis_functions.block_mass_matrix();
      FieldVector<double, 4> local_coeffs;
      FieldVector<double, 4> local_u;
      for (size_t jj = 0; jj < mass_matrix->num_rows / 4; ++jj) {
        for (size_t kk = 0; kk < 4; ++kk)
          local_u[kk] = u[4 * jj + kk];
        mass_matrix->block(jj).solve(local_coeffs, local_u);
        for (size_t kk = 0; kk < 4; ++kk)
          pn_coeffs[4 * jj + kk] = local_coeffs[kk];
      }
    }
  };

  template <bool anything>
  struct SolverHelper<false, anything>
  {
    static void solve(const MomentBasisType& basis_functions, const DynamicRangeType& u, DynamicRangeType& pn_coeffs)
    {
      const auto mass_matrix = basis_functions.mass_matrix();
      XT::LA::IstlRowMajorSparseMatrix<double> sparse_mass_matrix(mass_matrix, true);
      XT::LA::IstlDenseVector<double> pn_coeffs_istl(sparse_mass_matrix.rows());
      XT::LA::IstlDenseVector<double> u_istl(pn_coeffs_istl);
      for (size_t ii = 0; ii < u_istl.size(); ++ii)
        u_istl.set_entry(ii, u[ii]);
      XT::LA::solve(sparse_mass_matrix, u_istl, pn_coeffs_istl);
      for (size_t ii = 0; ii < pn_coeffs_istl.size(); ++ii)
        pn_coeffs[ii] = pn_coeffs_istl.get_entry(ii);
    }
  };

  static void run(const int quad_refinements, std::string testcase, std::string filename = "")
  {
    //******************* create grid and basisfunctions ***************************************
    auto grid_config = Dune::XT::Grid::cube_gridprovider_default_config();
    const auto grid_provider =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator());
    const auto grid_view = grid_provider.leaf_view();

    // ***************** choose testcase *********************
    RangeFieldType tau = 1e-12;
    RangeFieldType additional_quad_refinements = 0;
    RangeFieldType basis_quad_order = MomentBasisType::default_quad_order();
    RangeFieldType basis_quad_refinements = MomentBasisType::default_quad_refinements();
    std::function<RangeFieldType(const DomainType&, const bool)> psi;
    if (testcase == "Heaviside") {
      if (dimDomain != 1)
        DUNE_THROW(InvalidStateException,
                   "This is a 1-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      psi = [](const DomainType& v, const bool left) {
        return v[0] < 0 || (XT::Common::is_zero(v[0]) && left) ? RangeFieldType(5e-9) : RangeFieldType(1);
      };
    } else if (testcase == "Gauss1d") {
      if (dimDomain != 1)
        DUNE_THROW(InvalidStateException,
                   "This is a 1-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      const RangeFieldType sigma = 0.5;
      const RangeFieldType mu = 0.;
      psi = [sigma, mu](const DomainType& v, const bool) {
        return 1. / std::sqrt(2. * M_PI * std::pow(sigma, 2))
               * std::exp(std::pow(v[0] - mu, 2) / (-2 * std::pow(sigma, 2)));
      };
    } else if (testcase == "CrossingBeams1dSmooth") {
      if (dimDomain != 1)
        DUNE_THROW(InvalidStateException,
                   "This is a 1-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      const RangeFieldType factor = 1e5;
      const RangeFieldType norm = std::sqrt(M_PI) * std::erf(2 * std::sqrt(factor)) / (2 * std::sqrt(factor));
      psi = [factor, norm](const DomainType& v, const bool) {
        return std::max((std::exp(-factor * std::pow(v[0] - 1, 2)) + std::exp(-factor * std::pow(v[0] + 1, 2))) / norm,
                        5e-9);
      };
    } else if (testcase == "CrossingBeams1dDiscontinuous") {
      if (dimDomain != 1)
        DUNE_THROW(InvalidStateException,
                   "This is a 1-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      psi = [](const DomainType& v, const bool) {
        RangeFieldType height = 100;
        if (std::abs(v[0] + 1) < 1. / (2 * height) || std::abs(v[0] - 0.5) < 1. / (4 * height))
          return height;
        else
          return 5e-9;
      };
    } else if (testcase == "GaussOnSphere") {
      if (dimDomain != 3)
        DUNE_THROW(InvalidStateException,
                   "This is a 3-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      const RangeFieldType sigma = 0.5;
      const RangeFieldType mu = 1.;
      psi = [sigma, mu](const DomainType& v, const bool) {
        return 1. / (2. * M_PI * std::pow(sigma, 2) * (1. - std::exp(-2 / std::pow(sigma, 2))))
               * std::exp((std::pow(v[0] - mu, 2) + std::pow(v[1], 2) + std::pow(v[2], 2)) / (-2 * std::pow(sigma, 2)));
      };
    } else if (testcase == "SquareOnSphere") {
      if (dimDomain != 3)
        DUNE_THROW(InvalidStateException,
                   "This is a 3-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      psi = [](const DomainType& v, const bool) {
        if (v[0] > 0 && std::abs(v[1]) < 0.5 && std::abs(v[2]) < 0.5)
          return RangeFieldType(1);
        else
          return RangeFieldType(1e-8 / (4 * M_PI));
      };
    } else if (testcase == "CrossingBeams3d") {
      if (dimDomain != 3)
        DUNE_THROW(InvalidStateException,
                   "This is a 3-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      const RangeFieldType factor = 1e2;
      const RangeFieldType norm = M_PI / factor;
      const FieldVector<RangeFieldType, 3> center1{{1., 0., 0.}};
      const FieldVector<RangeFieldType, 3> center2{{0., 1., 0.}};
      psi = [factor, norm, center1, center2](const DomainType& v, const bool) {
        return std::max((std::exp(-factor * (v - center1).two_norm2()) + std::exp(-factor * (v - center2).two_norm2()))
                            / norm,
                        1e-8 / (4 * M_PI));
      };
    } else if (testcase == "CrossingBeams3d_2") {
      if (dimDomain != 3)
        DUNE_THROW(InvalidStateException,
                   "This is a 3-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      const RangeFieldType factor = 1e2;
      const RangeFieldType norm = M_PI / factor;
      const FieldVector<RangeFieldType, 3> center1{{std::sqrt(1. / 3.), std::sqrt(1. / 3.), std::sqrt(1. / 3.)}};
      const FieldVector<RangeFieldType, 3> center2{
          {std::sqrt(1. / (2 * M_PI)), -std::sqrt(0.5 - 1. / (4 * M_PI)), std::sqrt(0.5 - 1. / (4 * M_PI))}};
      psi = [factor, norm, center1, center2](const DomainType& v, const bool) {
        return std::max((std::exp(-factor * (v - center1).two_norm2()) + std::exp(-factor * (v - center2).two_norm2()))
                            / norm,
                        1e-8 / (4 * M_PI));
      };
    } else if (testcase == "CrossingBeams3dDiscontinuous") {
      if (dimDomain != 3)
        DUNE_THROW(InvalidStateException,
                   "This is a 3-dimensional test, but the basis functions are " + XT::Common::to_string(dimDomain)
                       + "-dimensional!");
      const FieldVector<RangeFieldType, 3> center1{{std::sqrt(1. / 3.), std::sqrt(1. / 3.), std::sqrt(1. / 3.)}};
      const FieldVector<RangeFieldType, 3> center2{
          {std::sqrt(1. / (2 * M_PI)), -std::sqrt(0.5 - 1. / (4 * M_PI)), std::sqrt(0.5 - 1. / (4 * M_PI))}};
      const RangeFieldType r = 0.01;
      const RangeFieldType r_squared = std::pow(r, 2);
      psi = [r_squared, center1, center2](const DomainType& v, const bool) {
        if ((v - center1).two_norm2() < r_squared || (v - center2).two_norm2() < r_squared)
          return 1 / (M_PI * r_squared);
        else
          return 1e-8 / (4 * M_PI);
      };
    } else {
      DUNE_THROW(NotImplemented, "Unknown testcase " + testcase + "!");
    }

    //******************* choose quadrature that is used for visualization and error calculation *******************
    const auto quadratures = QuadratureHelper<>::get(quad_refinements, additional_quad_refinements);
    std::shared_ptr<const MomentBasisType> basis_functions =
        std::make_shared<const MomentBasisType>(basis_quad_order, basis_quad_refinements);

    const auto u = basis_functions->get_moment_vector(psi);
    using MnFluxType = EntropyBasedFluxFunction<GridViewType, MomentBasisType>;
    MnFluxType mn_flux(grid_view, *basis_functions, tau, true);
    const auto mn_ret = mn_flux.get_alpha(u, true);
    if (mn_ret->second.second > 0.)
      std::cout << "Optimization problem regularized with r = " << mn_ret->second.second << " for the ansatz with "
                << dimRange << " moments!" << std::endl;
    const auto alpha = mn_ret->first;

    const auto entropy_string =
        "_"
        + std::string(MomentBasisType::entropy == EntropyType::MaxwellBoltzmann ? "MaxwellBoltzmann" : "BoseEinstein")
        + "_";
    const auto filename_mn = filename + entropy_string + basis_functions->mn_name();
    const auto filename_pn = filename + "_" + basis_functions->pn_name();
    std::ofstream mn_file(filename_mn + ".txt");
    std::ofstream pn_file(filename_pn + ".txt");
    const std::string mn_errors_filename = filename + entropy_string + basis_functions->short_id() + "m_errors.txt";
    const std::string pn_errors_filename = filename + basis_functions->short_id() + "p_errors.txt";
    std::ofstream mn_errors_file(mn_errors_filename, std::ios_base::app);
    std::ofstream pn_errors_file(pn_errors_filename, std::ios_base::app);
    if (is_empty(mn_errors_filename))
      mn_errors_file << "n l1error l2error linferror" << std::endl;
    if (is_empty(pn_errors_filename))
      pn_errors_file << "n l1error l2error linferror" << std::endl;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      mn_file << "v" << dd << " ";
      pn_file << "v" << dd << " ";
    }
    mn_file << "psi" << std::endl;
    pn_file << "psi" << std::endl;
    RangeFieldType l1norm(0), l2norm(0), linfnorm(0);
    RangeFieldType l1error_mn(0), l2error_mn(0), linferror_mn(0);
    RangeFieldType l1error_pn(0), l2error_pn(0), linferror_pn(0);
    DynamicRangeType pn_coeffs(dimRange);
    SolverHelper<>::solve(*basis_functions, u, pn_coeffs);
    const auto quadrature = XT::Data::merged_quadrature(quadratures);
    for (auto it = quadrature.begin(); it != quadrature.end(); ++it) {
      const auto& quad_point = *it;
      const auto v = quad_point.position();
      const auto weight = quad_point.weight();
      const auto basis = basis_functions->evaluate(v, it.first_index());
      const auto psi_mn_prod = std::inner_product(basis.begin(), basis.end(), alpha.begin(), 0.);
      const auto psi_mn = (MomentBasisType::entropy == EntropyType::MaxwellBoltzmann)
                              ? std::exp(psi_mn_prod)
                              : std::exp(psi_mn_prod) / (1 - std::exp(psi_mn_prod));
      const auto psi_pn = pn_coeffs * basis;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        mn_file << XT::Common::to_string(v[dd], 15) << " ";
        pn_file << XT::Common::to_string(v[dd], 15) << " ";
      }
      mn_file << XT::Common::to_string(psi_mn, 15) << std::endl;
      pn_file << XT::Common::to_string(psi_pn, 15) << std::endl;
      l1norm += std::abs(psi(v, basis_functions->is_negative(it))) * weight;
      l2norm += std::pow(psi(v, basis_functions->is_negative(it)), 2) * weight;
      linfnorm = std::max(std::abs(psi(v, basis_functions->is_negative(it))), linfnorm);
      l1error_mn += std::abs(psi_mn - psi(v, basis_functions->is_negative(it))) * weight;
      l2error_mn += std::pow(psi_mn - psi(v, basis_functions->is_negative(it)), 2) * weight;
      linferror_mn = std::max(std::abs(psi_mn - psi(v, basis_functions->is_negative(it))), linferror_mn);
      l1error_pn += std::abs(psi_pn - psi(v, basis_functions->is_negative(it))) * weight;
      l2error_pn += std::pow(psi_pn - psi(v, basis_functions->is_negative(it)), 2) * weight;
      linferror_pn = std::max(std::abs(psi_pn - psi(v, basis_functions->is_negative(it))), linferror_pn);
    }
    mn_errors_file << dimRange << " " << XT::Common::to_string(l1error_mn, 15) << " "
                   << XT::Common::to_string(l2error_mn, 15) << " " << XT::Common::to_string(linferror_mn, 15)
                   << std::endl;
    pn_errors_file << dimRange << " " << XT::Common::to_string(l1error_pn, 15) << " "
                   << XT::Common::to_string(l2error_pn, 15) << " " << XT::Common::to_string(linferror_pn, 15)
                   << std::endl;
  }
};

template <class MomentBasisType, class DiscreteFunctionType>
const size_t MomentApproximation<MomentBasisType, DiscreteFunctionType>::dimDomain;
template <class MomentBasisType, class DiscreteFunctionType>
const size_t MomentApproximation<MomentBasisType, DiscreteFunctionType>::dimRange;

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_HYPERBOLIC_MOMENT_APPROXIMATION_HH
