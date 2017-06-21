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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include <dune/gdt/local/fluxes/entropybased.hh>

#include "../base.hh"

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

template <class ImplementationType>
class KineticEquation : public ProblemBase<typename ImplementationType::EntityType,
                                           typename ImplementationType::DomainFieldType,
                                           ImplementationType::dimDomain,
                                           typename ImplementationType::StateType,
                                           typename ImplementationType::RangeFieldType,
                                           ImplementationType::dimRange>
{
protected:
  typedef ProblemBase<typename ImplementationType::EntityType,
                      typename ImplementationType::DomainFieldType,
                      ImplementationType::dimDomain,
                      typename ImplementationType::StateType,
                      typename ImplementationType::RangeFieldType,
                      ImplementationType::dimRange>
      BaseType;

public:
  using typename BaseType::RangeFieldType;

  KineticEquation(const ImplementationType& implementation)
    : BaseType(implementation.create_flux(),
               implementation.create_rhs(),
               implementation.create_initial_values(),
               implementation.create_boundary_values(),
               implementation.grid_config(),
               implementation.boundary_config(),
               implementation.CFL(),
               implementation.t_end())
  {
  }

  template <class BasisfunctionType, class GridLayerType>
  KineticEquation(const BasisfunctionType& basis_funcs, const GridLayerType grid_layer)
    : KineticEquation(ImplementationType(basis_funcs, grid_layer))
  {
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    return ImplementationType::default_grid_cfg();
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    return ImplementationType::default_grid_cfg();
  }

  bool has_non_zero_rhs() const
  {
    return true;
  }

  static std::string static_id()
  {
    return ImplementationType::static_id();
  }
}; // class KineticEquation<...>

template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          class RangeFieldImp,
          size_t rangeDim,
          size_t quadratureDim = domainDim>
class KineticTransportEquation
{
  typedef KineticTransportEquation ThisType;
  typedef XT::Functions::GlobalLambdaFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      GlobalLambdaFunctionType;
  typedef XT::Functions::GlobalLambdaFluxFunction<U_, 0, RangeFieldImp, rangeDim, 1> GlobalLambdaFluxFunctionType;

public:
  typedef BasisfunctionImp BasisfunctionType;
  typedef EntityImp EntityType;
  typedef DomainFieldImp DomainFieldType;
  typedef U_ StateType;
  typedef RangeFieldImp RangeFieldType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = rangeDim;

  typedef typename KineticEquation<ThisType>::FluxType FluxType;
  typedef typename KineticEquation<ThisType>::RhsType RhsType;
  typedef typename KineticEquation<ThisType>::InitialValueType InitialValueType;
  typedef typename KineticEquation<ThisType>::BoundaryValueType BoundaryValueType;

  typedef
      typename XT::Functions::AffineFluxFunction<EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange, 1>
          RhsAffineFunctionType;
  typedef typename XT::Functions::
      AffineFluxFunction<EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange, dimDomain>
          ActualFluxType;
  typedef XT::Functions::
      CheckerboardFunction<EntityImp, DomainFieldImp, dimDomain, RangeFieldImp, dimRange, 1, RhsAffineFunctionType>
          ActualRhsType;
  typedef XT::Functions::
      CheckerboardFunction<EntityImp, DomainFieldImp, dimDomain, RangeFieldImp, dimRange, 1, GlobalLambdaFunctionType>
          ActualInitialValueType;
  typedef GlobalLambdaFunctionType ActualBoundaryValueType;
  typedef typename RhsAffineFunctionType::FieldMatrixType MatrixType;
  typedef typename RhsAffineFunctionType::DomainType DomainType;
  typedef typename RhsAffineFunctionType::RangeType RangeType;
  typedef Dune::QuadratureRule<DomainFieldType, quadratureDim> QuadratureType;

  virtual ~KineticTransportEquation()
  {
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[1.0]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    XT::Common::Configuration boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  static QuadratureType default_quadrature(const XT::Common::Configuration& grid_cfg = default_grid_cfg())
  {
    std::vector<int> num_quad_cells = grid_cfg.get("num_quad_cells", std::vector<int>{2, 2, 2});
    size_t quad_order = grid_cfg.get("quad_order", 20);
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
                           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
                           const QuadratureType quadrature = QuadratureType(),
                           const RangeFieldType psi_vac = 5e-9)
    : basis_functions_(basis_functions)
    , grid_cfg_(grid_cfg)
    , boundary_cfg_(boundary_cfg)
    , psi_vac_(psi_vac)
    , quadrature_(quadrature)
  {
    if (quadrature_.empty())
      quadrature_ = default_quadrature(grid_cfg);
  }

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  virtual FluxType* create_flux() const
  {
    auto A = basis_functions_.mass_matrix_with_v();
    auto M_inv = basis_functions_.mass_matrix_inverse();
    for (size_t dd = 0; dd < dimDomain; ++dd)
      A[dd].rightmultiply(M_inv);
    return new ActualFluxType(A, typename ActualFluxType::RangeType(0));
  }

  virtual XT::Common::Parameter parameters() const
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{0}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.5}),
                                  std::make_pair("t_end", std::vector<double>{1.0}),
                                  std::make_pair("num_segments", std::vector<double>{1.})});
  }

  // RHS is (sigma_s/vol*G - sigma_t * I)u + Q<b>,
  // where sigma_t = sigma_s + sigma_a, G = <b><b>^T M^{-1} = <b>*c^T and
  // vol = <1> is the volume of the integration domain.
  virtual RhsType* create_rhs() const
  {
    const auto param = parameters();
    const auto num_segments = get_num_segments(param);
    auto sigma_a = param.get("sigma_a");
    auto sigma_s = param.get("sigma_s");
    auto Q = param.get("Q");
    const size_t num_regions =
        std::accumulate(num_segments.begin(), num_segments.end(), 1, [](auto a, auto b) { return a * b; });
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    auto sigma_t = sigma_a;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_t[ii] += sigma_s[ii];
    const RangeType basis_integrated = basis_functions_.integrated();
    const MatrixType M_inv = basis_functions_.mass_matrix_inverse();
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
    return new ActualRhsType(lower_left, upper_right, num_segments, affine_functions);
  } // ... create_rhs(...)

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is basis_integrated * psi_vac.
  virtual InitialValueType* create_initial_values() const
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    RangeType value = basis_functions_.integrated();
    const auto num_segments = get_num_segments(parameters());
    const size_t num_regions =
        std::accumulate(num_segments.begin(), num_segments.end(), 1, [](auto a, auto b) { return a * b; });
    value *= psi_vac_;
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    for (size_t ii = 0; ii < num_regions; ++ii)
      initial_vals.emplace_back([=](DomainType) { return value; }, 0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

  // Use a constant vacuum concentration basis_integrated * psi_vac as boundary value
  virtual BoundaryValueType* create_boundary_values() const
  {
    RangeType value = basis_functions_.integrated();
    value *= psi_vac_;
    return new ActualBoundaryValueType([=](const DomainType&) { return value; }, 0);
  } // ... create_boundary_values()

  virtual RangeFieldType CFL() const
  {
    return parameters().get("CFL")[0];
  }

  virtual RangeFieldType t_end() const
  {
    return parameters().get("t_end")[0];
  }

  virtual XT::Common::Configuration grid_config() const
  {
    return grid_cfg_;
  }

  virtual XT::Common::Configuration boundary_config() const
  {
    return boundary_cfg_;
  }

  static std::string static_id()
  {
    return "kinetictransportequation";
  }

  const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

protected:
  FieldVector<size_t, dimDomain> get_num_segments(const XT::Common::Parameter& param) const
  {
    FieldVector<size_t, dimDomain> num_segments;
    std::vector<double> stdvec_num_segments = param.get("num_segments");
    for (size_t dd = 0; dd < dimDomain; ++dd)
      num_segments[dd] = size_t(stdvec_num_segments[dd] + 0.1); // 0.1 ensures correct conversion in all cases
    return num_segments;
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

  const BasisfunctionType& basis_functions_;
  const XT::Common::Configuration grid_cfg_;
  const XT::Common::Configuration boundary_cfg_;
  const RangeFieldType psi_vac_;
  QuadratureType quadrature_;
}; // class KineticTransportEquation<...>


template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          class RangeFieldImp,
          size_t rangeDim>
class KineticFokkerPlanckEquation : public KineticTransportEquation<BasisfunctionImp,
                                                                    EntityImp,
                                                                    DomainFieldImp,
                                                                    domainDim,
                                                                    U_,
                                                                    RangeFieldImp,
                                                                    rangeDim>
{
  typedef KineticTransportEquation<BasisfunctionImp, EntityImp, DomainFieldImp, domainDim, U_, RangeFieldImp, rangeDim>
      BaseType;

public:
  using typename BaseType::BasisfunctionType;
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

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  KineticFokkerPlanckEquation(const BasisfunctionType& basis_functions,
                              const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                              const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
                              const RangeFieldType psi_vac = 1e-4)
    : BaseType(basis_functions, grid_cfg, boundary_cfg, psi_vac)
  {
  }


  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{4}),
                                  std::make_pair("T", std::vector<double>{0}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.5}),
                                  std::make_pair("t_end", std::vector<double>{4.0}),
                                  std::make_pair("num_segments", std::vector<double>{1.})});
  }

  // RHS is (-\sigma_a*I + 0.5*T*S M^{-1}) u + Q<b>
  virtual RhsType* create_rhs() const override
  {
    const auto param = parameters();
    const auto num_segments = get_num_segments(param);
    const std::vector<RangeFieldType> sigma_a = param.get("sigma_a");
    const std::vector<RangeFieldType> T = param.get("T");
    const std::vector<RangeFieldType> Q = param.get("Q");
    const size_t num_regions =
        std::accumulate(num_segments.begin(), num_segments.end(), 1, [](auto a, auto b) { return a * b; });
    assert(sigma_a.size() == T.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const RangeType basis_integrated = basis_functions_.integrated();
    const MatrixType M_inv = basis_functions_.mass_matrix_inverse();
    const MatrixType S = basis_functions_.S();
    MatrixType I(0);
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
    return new ActualRhsType(lower_left, upper_right, num_segments, affine_functions);
  } // ... create_rhs(...)

protected:
  using BaseType::get_num_segments;
  using BaseType::basis_functions_;
  using BaseType::grid_cfg_;
}; // class KineticFokkerPlanckEquation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH
