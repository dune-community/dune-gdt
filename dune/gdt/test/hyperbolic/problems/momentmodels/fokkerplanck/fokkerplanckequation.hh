// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCKEQUATION_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCKEQUATION_HH

#include <dune/xt/common/parameter.hh>

#include "../kinetictransport/kinetictransportequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


/** Testcase for the \f$P_n\f$ moment approximation of the Fokker-Planck equation in one dimension
 * \f[
 * \partial_t \psi(t,x,v) + v * \partial_x \psi(t,x,v) + \sigma_a(x)*\psi(t,x,v) = 0.5*T(x)*\Delta_v \psi(t,x,v) + Q(x),
 * \f]
 * where \f$\psi: [0,T] \times [x_l, x_dimRange] \times [-1, 1] \to \mathbb{R}\f$,
 * \f$Delta_v \psi = \partial_v( (1-v^2)\partial_v \psi)\f$ is the Laplace-Beltrami operator and
 * \f$\sigma_a, T, Q: [x_l, x_dimRange] \to \mathbb{R}\f$ are the absorption coefficient, the transport coefficient and
 * a
 * source, respectively.
 * The \f$P_n\f$ model approximates the solution of the Fokker-Planck equation by an ansatz
 * \f[
 * \psi(t,x,v) = \sum \limits_{l=0}^n u_i(t,x)\phi_i(v)
 * \f]
 * where the \f$\phi_i, \, i=0,\ldots,n$ are suitable basis functions of (a subset of) the function space on
 * \f$[-1, 1]\f$ that \f$\psi\f$ lives in. n is called the moment order. Usually, the \f$\phi_i\f$ are chosen as the
 * Legendre polynomials up to order n.
 * Once suitable basis functions are found, a Galerkin semidiscretization in v is done, so the \f$\phi_i\f$ are also
 * taken as basis for the test space. This results in an equation of the form
 * \f[
 * M \partial_t u + B \partial_x u = q - (\sigma_a*M + 0.5*T*S) u,
 * \f]
 *  where \f$u = (u_1, \ldots, u_n)^T\f$, \f$M, B, S \in \mathbb{R}^{n\times n}\f$ with
 * \f$ M_{ji} = (\phi_i, \phi_j)_v\f$, \f$B{ji} = (v*\phi_i, \phi_j)_v\f$,
 * \f$S_{ji} = ((1-v^2)\partial_v \phi_i, \partial_v \phi_j)_v\f$ and
 * \f$q_i(x) = (Q(x), \phi_i(v))_v = Q(x) (1, \phi_i(v))_v\f$. Here, \f$(a,b)_v = \int \limits_{-1}^1 a(v)b(v)dv\f$
 * is the \f$L^2\f$ inner product with respect to \f$v\f$.
 * In the following, we rescale \f$u\f$ s.t. \f$(\psi(t,x,v),\phi_i(v))_v = u_i(t,x)\f$ if the ansatz holds. Provided
 * the \f$ \phi_i \f$ are an orthogonal basis, \f$M\f$ is invertible and the rescaling corresponds to a multiplication
 * of \f$u\f$ by \f$M^{-1}\f$ from the left, giving the equation
 * \f[
 * \partial_t u + B M^{-1} \partial_x u = q - (\sigma_a*I_{n\times n} + 0.5*T*S M^{-1}) u.
 * \f]
 * For details on the parameters of the test cases implemented here (OneBeam, TwoBeams, TwoPulses, RectangularIC,
 * SourceBeam) see Schneider, Alldredge, Frank, Klar, "Higher Order Mixed-Moment Approximations for the
 * Fokker-Planck Equation in One Space Dimension", SIAM J. Appl. Math., 74(4), 1087–1114
 */
template <class BasisfunctionImp, class GridLayerImp, class U_>
class FokkerPlanckEquation : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>
{
  typedef KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_> BaseType;

public:
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using BaseType::dimDomain;
  static_assert(BaseType::dimDomain == 1, "Not implemented for dimDomain > 1!");
  using BaseType::dimRange;
  using typename BaseType::RhsType;
  using typename BaseType::ActualRhsType;
  using typename BaseType::RhsAffineFunctionType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  FokkerPlanckEquation(const BasisfunctionType& basis_functions,
                       const GridLayerType grid_layer,
                       const QuadratureType quadrature = QuadratureType(),
                       const size_t num_segments = 1,
                       const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                       const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
                       const RangeFieldType psi_vac = 1e-4)
    : BaseType(basis_functions,
               grid_layer,
               quadrature,
               {num_segments},
               grid_cfg,
               boundary_cfg,
               psi_vac,
               XT::Common::ParameterType({std::make_pair("sigma_a", get_num_regions({num_segments})),
                                          std::make_pair("T", get_num_regions({num_segments})),
                                          std::make_pair("Q", get_num_regions({num_segments})),
                                          std::make_pair("CFL", 1),
                                          std::make_pair("t_end", 1)}))
  {
  }

  using BaseType::parse_parameter;
  using BaseType::parameters;

  // RHS is (-\sigma_a*I + 0.5*T*S M^{-1}) u + Q<b>
  virtual RhsType* create_rhs() const override
  {
    const auto param = parse_parameter(parameters());
    const std::vector<RangeFieldType> sigma_a = param.get("sigma_a");
    const std::vector<RangeFieldType> T = param.get("T");
    const std::vector<RangeFieldType> Q = param.get("Q");
    const size_t num_regions = get_num_regions(num_segments_);
    assert(sigma_a.size() == T.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const RangeType basis_integrated = basis_functions_.integrated();
    const MatrixType M_inv = basis_functions_.mass_matrix_inverse();
    const MatrixType S = basis_functions_.S();
    MatrixType I(dimRange, dimRange, 0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType K = S;
    K.rightmultiply(M_inv);

    std::vector<RhsAffineFunctionType> affine_functions;
    for (size_t ii = 0; ii < num_regions; ++ii) {
      MatrixType K_scaled = K;
      K_scaled *= T[ii] / 2.;
      MatrixType I_scaled = I;
      I_scaled *= sigma_a[ii];
      MatrixType A = K_scaled;
      A -= I_scaled;
      RangeType b = basis_integrated;
      b *= Q[ii];
      affine_functions.emplace_back(A, b);
    } // ii
    return new ActualRhsType(lower_left, upper_right, num_segments_, affine_functions);
  } // ... create_rhs(...)

protected:
  using BaseType::num_segments_;
  using BaseType::get_num_regions;
  using BaseType::basis_functions_;
  using BaseType::grid_cfg_;
}; // class FokkerPlanckEquation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCKEQUATION_HH
