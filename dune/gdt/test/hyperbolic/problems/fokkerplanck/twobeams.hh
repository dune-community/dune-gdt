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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/global.hh>

#include "../base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


// forwards
template <class BasisFunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class TwoBeamsPn;


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
template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          class RangeFieldImp,
          size_t rangeDim>
class KineticProblemTraitsBase
{
  typedef XT::Functions::GlobalLambdaFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      GlobalLambdaFunctionType;
  typedef XT::Functions::
      GlobalLambdaMetaFunction<EntityImp, DomainFieldImp, domainDim, U_, 0, RangeFieldImp, rangeDim, 1>
          GlobalLambdaMetaFunctionType;

public:
  typedef BasisfunctionImp BasisFunctionType;
  typedef EntityImp EntityType;
  typedef DomainFieldImp DomainFieldType;
  typedef U_ StateType;
  typedef RangeFieldImp RangeFieldType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = rangeDim;

  typedef
      typename XT::Functions::AffineMetaFunction<EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange, 1>
          RhsAffineFunctionType;
  typedef typename XT::Functions::
      AffineMetaFunction<EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange, dimDomain>
          ActualFluxType;
  typedef XT::Functions::
      CheckerboardFunction<EntityImp, DomainFieldImp, dimDomain, RangeFieldImp, dimRange, 1, RhsAffineFunctionType>
          ActualRhsType;
  typedef XT::Functions::
      CheckerboardFunction<EntityImp, DomainFieldImp, dimDomain, RangeFieldImp, dimRange, 1, GlobalLambdaFunctionType>
          ActualInitialValueType;
  typedef GlobalLambdaFunctionType ActualBoundaryValueType;
}; // class KineticProblemTraitsBase<...>

template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          class RangeFieldImp,
          size_t rangeDim>
class TwoBeamsPnTraits : public KineticProblemTraitsBase<BasisfunctionImp,
                                                         EntityImp,
                                                         DomainFieldImp,
                                                         domainDim,
                                                         U_,
                                                         RangeFieldImp,
                                                         rangeDim>
{
  typedef TwoBeamsPn<BasisfunctionImp, EntityImp, DomainFieldImp, domainDim, U_, RangeFieldImp, rangeDim> derived_type;
}; // class TwoBeamPnTraits<...>

template <class Traits>
class KineticTransportEquation : public ProblemBase<typename Traits::EntityType,
                                                    typename Traits::DomainFieldType,
                                                    Traits::dimDomain,
                                                    typename Traits::StateType,
                                                    0,
                                                    typename Traits::RangeFieldType,
                                                    Traits::dimRange>
{
  typedef ProblemBase<typename Traits::EntityType,
                      typename Traits::DomainFieldType,
                      dimDomain,
                      typename Traits::StateType,
                      0,
                      typename Traits::RangeFieldType,
                      dimRange>
      BaseType;
  typedef Traits::derived_type ProblemType;

public:
  using typename BaseType::BasisfunctionType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  typedef typename Traits::ActualFluxType ActualFluxType;
  typedef typename Traits::ActualRhsType ActualRhsType;
  typedef typename Traits::ActualInitialValueType ActualInitialValueType;
  typedef typename Traits::ActualBoundaryValueType ActualBoundaryValueType;

  static XT::Common::Configuration default_grid_config()
  {
    return ProblemType::default_grid_config();
  }

  static XT::Common::Configuration default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  template <class... Args>
  KineticTransportEquation(const BasisfunctionType& basis_functions,
                           const XT::Common::Configuration& grd_cfg = ProblemType::default_grid_cfg(),
                           const XT::Common::Configuration& bnd_cfg = ProblemType::default_boundary_info_cfg(),
                           Args&&... args)
    : BaseType(ProblemType::create_flux(basis_functions, grid_cfg, std::forward<Args>(args)...),
               ProblemType::create_rhs(basis_functions, grid_cfg, std::forward<Args>(args)...),
               ProblemType::create_initial_values(basis_functions, grid_cfg, std::forward<Args>(args)...),
               ProblemType::create_boundary_values(basis_functions, grid_cfg, std::forward<Args>(args)...),
               grd_cfg,
               bnd_cfg)
  {
  }

protected:
  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  static FluxType* create_flux(const BasisFunctionType& basis_functions, const XT::Common::Configuration& /*grd_cfg*/)
  {
    auto A = basis_functions.mass_matrix_with_v();
    auto M_inv = basis_functions.mass_matrix_inverse();
    for (size_t dd = 0; dd < dimDomain; ++dd)
      A[dd].rightmultiply(M_inv);
    return new ActualFluxType(A, RangeType(0));
  }

  static RangeFieldType volume()
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

  // RHS is (sigma_s/vol*G - sigma_t * I)u + Q<b>,
  // where sigma_t = sigma_s + sigma_a, G = <b><b>^T M^{-1} = <b>*c^T and
  // vol = <1> is the volume of the integration domain.
  static FluxType* create_rhs(const BasisfunctionType& basis_functions, const XT::Common::Configuration& grd_cfg)
  {
    typedef typename ProblemTraits::RhsAffineFunctionType AffineFunctionType;
    typedef typename AffineFunctionType::FieldMatrixType MatrixType;
    const FieldVector<size_t, 3> num_elements = ProblemType::num_elements();
    const std::vector<RangeFieldType> sigma_a = ProblemType::sigma_a();
    const std::vector<RangeFieldType> sigma_s = ProblemType::sigma_s();
    const std::vector<RangeFieldType> Q = TestcaseType::Q();
    const size_t num_regions = std::accumulate(num_elements.begin(), num_elements.end());
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grd_cfg["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grd_cfg["upper_right"]);
    const auto sigma_t = sigma_a;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_t[ii] += sigma_s[ii];
    const RangeType basis_integrated = basis_functions.integrated();
    const MatrixType M_inv = basis_functions.mass_matrix_inverse();
    RangeType c(0);
    M_inv.mtv(basis_integrated, c);
    MatrixType I(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType G(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        G[rr][cc] = basis_integrated[rr] * c[cc];
    const auto vol = volume();

    std::vector<typename ProblemTraits::RhsAffineFunctionType> affine_functions;
    for (size_t ii = 0; ii < num_regions; ++ii) {
      MatrixType G_scaled = G;
      G_scaled *= sigma_s[ii] / vol;
      MatrixType I_scaled = I;
      I_scaled *= sigma_t[ii];
      MatrixType A = G_scaled;
      A -= I_scaled;
      RangeType b = basis_integrated;
      b *= Q[ii];
      affine_functions.emplace_back(A, b);
    } // ii
    return std::make_unique<ActualRhsType>(lower_left, upper_right, num_elements, affine_functions);
  } // ... create_rhs(...)
}; // class KineticTransportEquation<...>

template <class Traits>
class KineticFokkerPlanckEquation : public KineticTransportEquation<Traits>
{
  typedef KineticTransportEquation<Traits> BaseType;
  using typename BaseType::ProblemType;

public:
  using typename BaseType::BasisfunctionType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  static_assert(BaseType::dimDomain == 1, "Not implemented for dimDomain > 1!");
  static const size_t dimRange = BaseType::dimRange;
  using typename BaseType::RhsType;
  typedef typename Traits::ActualRhsType ActualRhsType;
  using typename BaseType::BasisFunctionType;

  template <class... Args>
  KineticFokkerPlanckEquation(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

protected:
  // RHS is (-\sigma_a*I + 0.5*T*S M^{-1}) u + Q<b>
  static RhsType* create_rhs(const BasisfunctionType& basis_functions, const XT::Common::Configuration& grd_cfg)
  {
    typedef typename Traits::RhsAffineFunctionType AffineFunctionType;
    typedef typename AffineFunctionType::FieldMatrixType MatrixType;
    const FieldVector<size_t, 3> num_elements = ProblemType::num_elements();
    const std::vector<RangeFieldType> sigma_a = ProblemType::sigma_a();
    const std::vector<RangeFieldType> T = ProblemType::T();
    const std::vector<RangeFieldType> Q = TestcaseType::Q();
    const size_t num_regions = std::accumulate(num_elements.begin(), num_elements.end());
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grd_cfg["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grd_cfg["upper_right"]);
    const RangeType basis_integrated = basis_functions.integrated();
    const MatrixType M_inv = basis_functions.mass_matrix_inverse();
    const S = basis_functions.S();
    MatrixType I(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType K = S;
    K.rightmultiply(M_inv);

    std::vector<AffineFunctionType> affine_functions;
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
    return std::make_unique<ActualRhsType>(lower_left, upper_right, num_elements, affine_functions);
  } // ... create_rhs(...)
}; // class KineticFokkerPlanckEquation<...>


template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange,
          class Traits =
              TwoBeamsPnTraits<BasisfunctionImp, EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange>>
class TwoBeamsPn : public KineticFokkerPlanckEquation<Traits>
{
  typedef typename KineticTransportEquation<Traits> BaseType;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;

  template <class... Args>
  TwoBeamsPn(Args... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static FieldVector<size_t, dimDomain> num_elements()
  {
    return FieldVector<size_t, dimDomain>(1);
  }

  static RangeFieldType sigma_a()
  {
    return 4;
  }

  static RangeFieldType T()
  {
    return 0;
  }

  static RangeFieldType Q()
  {
    return 0;
  }

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is 2 psi_vac if n == 0 and 0 else.
  static InitialValueType* create_initial_values(const BasisfunctionImp& basis_functions,
                                                 const XT::Common::Configuration& grd_cfg,
                                                 const RangeFieldType psi_vac = 1e-4)
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grd_cfg["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grd_cfg["upper_right"]);
    RangeType value(0);
    value[0] = 2 * psi_vac;
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back([=](const DomainType&) { return value; });
    return new ActualInitialValueType(
        lower_left, upper_right, TwoBeamsPn::num_elements(), initial_vals, "initial_values");
  } // ... create_initial_values()

  // boundary value of kinetic equation is 100*delta(v-1) at x = 0 and 100*delta(v+1) at x = 1,
  // so k-th component of boundary value has to be 50*\phi_k(1) at x = 0 and 50*\phi_k(-1) at x = 1.
  // Model with function(x) = 50*\phi_k(-1)*x + 50*\phi_k(1)*(1-x).
  static BoundaryValueType* create_boundary_values(const BasisfunctionImp& basis_functions,
                                                   const XT::Common::Configuration& /*grd_cfg*/)
  {
    const auto basis_evaluated_at_one = basis_functions.evaluate(DomainType(1));
    const auto basis_evaluated_at_minus_one = basis_functions.evaluate(DomainType(-1));
    return new ActualBoundaryValueType(
        [=](const DomainType& x) {
          RangeType ret = basis_evaluated_at_minus_one;
          ret *= x[0] * 50;
          RangeType summand2 = basis_evaluated_at_one;
          summand2 *= (1 - x[0]) * 50;
          ret += summand_2;
          return ret;
        },
        1);
  } // ... create_boundary_value_config()
};


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
