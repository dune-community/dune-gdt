// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container.hh>

#include "../default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** Testcase for the \f$P_n\f$ moment approximation of the Fokker-Planck equation in one dimension
 * \f[
 * \partial_t \psi(t,x,v) + v * \partial_x \psi(t,x,v) + \sigma_a(x)*\psi(t,x,v) = 0.5*T(x)*\Delta_v \psi(t,x,v) + Q(x),
 * \f]
 * where \f$\psi: [0,T] \times [x_l, x_r] \times [-1, 1] \to \mathbb{R}\f$,
 * \f$Delta_v \psi = \partial_v( (1-v^2)\partial_v \psi)\f$ is the Laplace-Beltrami operator and
 * \f$\sigma_a, T, Q: [x_l, x_r] \to \mathbb{R}\f$ are the absorption coefficient, the transport coefficient and a
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
 * M \partial_t u + D \partial_x u = q - (\sigma_a*M + 0.5*T*S) u,
 * \f]
 *  where \f$u = (u_1, \ldots, u_n)^T\f$, \f$M, D, S \in \mathbb{R}^{n\times n}\f$ with
 * \f$ M_{ji} = (\phi_i, \phi_j)_v\f$, \f$D_{ji} = (v*\phi_i, \phi_j)_v\f$,
 * \f$S_{ji} = ((1-v^2)\partial_v \phi_i, \partial_v \phi_j)_v\f$ and
 * \f$q_i(x) = (Q(x), \phi_i(v))_v = Q(x) (1, \phi_i(v))_v\f$. Here, \f$(a,b)_v = \int \limits_{-1}^1 a(v)b(v)dv\f$
 * is the \f$L^2\f$ inner product with respect to \f$v\f$.
 * In the following, we rescale \f$u\f$ s.t. \f$(\psi(t,x,v),\phi_i(v))_v = u_i(t,x)\f$ if the ansatz holds. Provided
 * the \f$ \phi_i \f$ are an orthogonal basis, \f$M\f$ is invertible and the rescaling corresponds to a multiplication
 * of \f$u\f$ by \f$M^{-1}\f$ from the left, giving the equation
 * \f[
 * \partial_t u + D M^{-1} \partial_x u = q - (\sigma_a*I_{n\times n} + 0.5*T*S M^{-1}) u.
 * \f]
 * For details on the parameters of the test cases implemented here (OneBeam, TwoBeams, TwoPulses, RectangularIC,
 * SourceBeam) see Schneider, Alldredge, Frank, Klar, "Higher Order Mixed-Moment Approximations for the
 * Fokker-Planck Equation in One Space Dimension", SIAM J. Appl. Math., 74(4), 1087–1114
 */
template <class E, class D, size_t d, class R, size_t momentOrder>
class TwoBeams : public Default<E, D, d, R, momentOrder + 1>
{
  typedef TwoBeams<E, D, d, R, momentOrder> ThisType;
  typedef Default<E, D, d, R, momentOrder + 1> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DummyEntityType;
  typedef typename Dune::Stuff::Functions::Affine<DummyEntityType, R, dimRange, R, dimRange, dimDomain>
      FluxAffineFunctionType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux<FluxAffineFunctionType, E, D, d, R, dimRange, 1>
      DefaultFluxType;
  typedef typename DefaultFluxType::RangeType RangeType;
  typedef typename DefaultFluxType::FluxRangeType FluxRangeType;
  typedef typename FluxAffineFunctionType::FieldMatrixType MatrixType;
  using typename BaseType::DefaultInitialValueType;
  typedef typename DS::Functions::Affine<DummyEntityType, R, dimRange, R, dimRange, 1> RHSAffineFunctionType;
  typedef typename DS::Functions::FunctionCheckerboard<RHSAffineFunctionType, E, D, d, R, dimRange, 1>
      RHSCheckerboardFunctionType;
  typedef typename Dune::GDT::CheckerboardBasedRhsEvaluationFluxInterface<RHSCheckerboardFunctionType, E, D, d, R,
                                                                          dimRange, 1> DefaultRHSType;
  typedef typename DefaultRHSType::DomainType DomainType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".twobeams";
  }

  std::string type() const override
  {
    return BaseType::type() + ".twobeams";
  }


  static std::string short_id()
  {
    return "2Beams";
  }

protected:
  class GetData
  {
  public:
    //    typedef D                                                            VelocityFieldImp;
    //    typedef typename Dune::YaspGrid< dimDomain, Dune::EquidistantOffsetCoordinates< double, dimDomain > >
    //    VelocityGridType;
    //    typedef Dune::Stuff::Grid::Providers::Cube< VelocityGridType >                    VelocityGridProviderType;
    //    typedef typename VelocityGridType::LeafGridView                                   VelocityGridViewType;
    //    typedef typename VelocityGridType::template Codim< 0 >::Entity                    VelocityEntityType;
    //    typedef typename DS::LocalizableFunctionInterface< VelocityEntityType,
    //                                                       VelocityFieldImp, dimDomain,
    //                                                       R, 1, 1 >          VelocityFunctionType;
    //    typedef typename DS::Functions::Expression< VelocityEntityType,
    //                                                VelocityFieldImp, dimDomain,
    //                                                R, 1, 1 >                 VelocityExpressionFunctionType;

    typedef typename Dune::Stuff::LA::CommonDenseVector<R> VectorType;
    //    typedef typename Dune::GDT::CgSpaceProvider< VelocityGridType,
    //                                                    DSG::ChooseLayer::leaf,
    //                                                    Dune::GDT::ChooseSpaceBackend::pdelab,
    //                                                    1, R, 1, 1 >          CGProviderType;
    //    typedef typename CGProviderType::Type                                             CGSpaceType;
    //    typedef Dune::GDT::DiscreteFunction< CGSpaceType, VectorType >                    CGFunctionType;
    //    typedef typename DS::Functions::Checkerboard< typename VelocityGridType::template Codim< 0 >::Entity,
    //                                                  D, dimDomain,
    //                                                  R, 1, 1 >               CGJacobianType;
    static const int precision = 15; // precision for to_string

    static void set_filename(const std::string filename)
    {
      filename_ = filename;
    }

    static void set_exact_legendre(const bool use_exact_legendre)
    {
      exact_legendre_ = use_exact_legendre;
    }

    static const bool& exact_legendre()
    {
      return exact_legendre_;
    }

    //    static const std::vector< CGFunctionType >& basefunctions()
    //    {
    //      calculate();
    //      return basefunctions_;
    //    }

    //    static std::shared_ptr< const VelocityGridViewType > velocity_grid_view()
    //    {
    //      calculate();
    //      return velocity_grid_view_;
    //    }

    // q - (sigma_a + T/2*S*M^(-1))*u = Q(x)*base_integrated() - (sigma_a(x)*I_{nxn} + T(x)/2*S*M_inverse)*u = q - A*u
    // here, T = 0, Q = 0, sigma_a = 4
    static void create_rhs_values(ConfigType& rhs_config)
    {
      std::string A_str = "[";
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (rr > 0)
          A_str += "; ";
        for (size_t cc = 0; cc < dimRange; ++cc) {
          if (cc > 0)
            A_str += " ";
          if (rr == cc)
            A_str += "-4";
          else
            A_str += "0";
        }
      }
      A_str += "]";
      rhs_config["A.0"] = A_str;
      rhs_config["b.0"] = DSC::to_string(FluxRangeType(0));
    } // ... create_rhs_values(...)

    // flux matrix is D*M^(-1)
    // for legendre polynomials, using a recursion relation gives D*M^(-1)[cc][rr] = rr/(2*rr + 1)       if cc == rr - 1
    //                                                                             = (rr + 1)/(2*rr + 1) if cc == rr + 1
    //                                                                             = 0                   else
    static std::string create_flux_matrix()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += "; ";
          for (size_t cc = 0; cc < dimRange; ++cc) {
            if (cc > 0)
              str += " ";
            if (cc == rr - 1)
              str += DSC::to_string(double(rr) / (2.0 * double(rr) + 1.0), precision);
            else if (cc == rr + 1)
              str += DSC::to_string((double(rr) + 1.0) / (2.0 * double(rr) + 1.0), precision);
            else
              str += "0";
          }
        }
        str += "]";
        return str;
      } else {
        MatrixType D_M_inverse(M_inverse());
        return DSC::to_string(D_M_inverse.leftmultiply(DD()), precision);
      }
    } // ... create_flux_matrix()

    // initial value of kinetic equation is constant 10^(-4), thus initial value of the k-th component of the
    // moment vector is 10^(-4)*base_integrated_k.
    // For Legendre polynomials, this is 0.0002 if k == 0 and 0 else
    static std::string create_initial_values()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          if (rr == 0)
            str += "0.0002";
          else
            str += "0";
        }
        str += "]";
        return str;
      } else {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += DSC::to_string(0.0001 * base_integrated()[rr], precision);
        }
        str += "]";
        return str;
      }
    } // ... create_initial_values()

    // boundary value of kinetic equation is 100*delta(v-1) at x = 0 and 100*delta(v+1) at x = 1,
    // so k-th component of boundary value has to be 50*\phi_k(1) at x = 0 and 50*\phi_k(-1) at x = 1
    // simulate with function(x) = 50*((\phi_k(-1) - \phi_k(1))*x + \phi_k(1))
    // For Legendre polynomials, this is [50 50 50 ...] at x = 0 and [50 -50 50 -50 ... ] at x = 1
    // simulate with function(x) = 50*((-1)^n - 1)*x + 1)
    static std::string create_boundary_values()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += "50*(" + DSC::to_string(((1.0 - 2.0 * (rr % 2)) - 1.0), precision) + "*x[0]+1)";
        }
        str += "]";
        return str;
      } else {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += DSC::to_string(50 * (basefunctions_values_at_minusone()[rr] - basefunctions_values_at_plusone()[rr]),
                                precision)
                 + "*x[0]+" + DSC::to_string(50 * basefunctions_values_at_plusone()[rr], precision);
        }
        str += "]";
        return str;
      }
    } // ... create_boundary_values()

  protected:
    static const MatrixType& M()
    {
      calculate();
      return M_;
    }

    static const MatrixType& DD()
    {
      calculate();
      return D_;
    }

    static const MatrixType& S()
    {
      calculate();
      return S_;
    }

    static const MatrixType& M_inverse()
    {
      calculate();
      return M_inverse_;
    }

    static const RangeType& base_integrated()
    {
      calculate();
      return base_integrated_;
    }

    static const RangeType& basefunctions_values_at_minusone()
    {
      calculate();
      return basefunctions_values_at_minusone_;
    }

    static const RangeType& basefunctions_values_at_plusone()
    {
      calculate();
      return basefunctions_values_at_plusone_;
    }

    static const RangeType& onebeam_left_boundary_values()
    {
      calculate();
      return onebeam_left_boundary_values_;
    }

  private:
    static void calculate()
    {
      //      if (!is_calculated_) {
      //        // get basis functions
      //        std::ifstream basefunction_file(filename_);
      //        // get grid points from first line
      //        std::string grid_points;
      //        getline(basefunction_file, grid_points);
      //        // get values of basefunctions at the DoFs, each line is for one basis function
      //        std::vector< std::vector< std::string > > basefunction_values(dimRange);
      //        for (size_t ii = 0; ii < dimRange; ++ii) {
      //          std::string line;
      //          std::getline(basefunction_file, line);
      //          basefunction_values[ii] = DSC::tokenize(line,
      //                                                  ",",
      //                                                  boost::algorithm::token_compress_mode_type::token_compress_on);
      //        }

      //        //create grid for the velocity space, the size has to be the number of DOFs - 1
      //        ConfigType velocity_grid_config;
      //        velocity_grid_config["type"] = "provider.cube";
      //        velocity_grid_config["lower_left"] = "[-1.0]";
      //        velocity_grid_config["upper_right"] = "[1.0]";
      //        velocity_grid_config["num_elements"] = "[" + DSC::to_string(basefunction_values[0].size() - 1) + "]";
      //        VelocityGridProviderType velocity_grid_provider =
      //        *(VelocityGridProviderType::create(velocity_grid_config));
      //        velocity_grid_ = velocity_grid_provider.grid_ptr();

      //        // make CG Space with polOrder 1 and DiscreteFunctions for the basis functions
      //        velocity_grid_view_ = std::make_shared< VelocityGridViewType >(velocity_grid_->leafGridView());

      //        CGSpaceType cg_space = CGProviderType::create(velocity_grid_provider);

      //        // create basefunctions from values at the DOFs
      //        for (size_t ii = 0; ii < dimRange; ++ii) {
      //          VectorType basefunction_ii_values(velocity_grid_view_->size(0) + 1);
      //          for (size_t jj = 0; jj < basefunction_values[ii].size(); ++jj) {
      //            basefunction_ii_values[jj] = DSC::from_string< R >(basefunction_values[ii][jj]);
      //          }
      //          basefunctions_values_at_minusone_[ii] = basefunction_ii_values[0];
      //          basefunctions_values_at_plusone_[ii] = basefunction_ii_values[velocity_grid_view_->size(0)];
      //          basefunctions_.emplace_back(CGFunctionType(cg_space,
      //                                                    basefunction_ii_values,
      //                                                    "Basefunction " + DSC::to_string(ii)));
      //        }

      //        // get jacobians of basefunctions. jacobians are piecewise constant, so use Checkerboard as
      //        CGJacobianType
      //        std::vector< CGJacobianType > basefunction_jacobians;
      //        std::vector< std::vector< typename CGJacobianType::RangeType > >
      //        basefunction_jacobians_values(dimRange);
      //        for (size_t ii = 0; ii < dimRange; ++ii) {
      //          basefunction_jacobians_values[ii].resize(velocity_grid_view_->size(0));
      //        }
      //        const auto it_end = velocity_grid_view_->template end< 0 >();
      //        for (auto it = velocity_grid_view_->template begin< 0 >(); it != it_end; ++it) {
      //          const auto& entity = *it;
      //          for (size_t ii = 0; ii < dimRange; ++ii) {
      //            // basefunctions_[ii].jacobian(..) gives a 1x1 FieldMatrix
      //            basefunction_jacobians_values[ii][velocity_grid_view_->indexSet().index(entity)]
      //                =
      //                (basefunctions_[ii].local_function(entity)->jacobian(entity.geometry().local(entity.geometry().center())))[0][0];
      //          }
      //        }

      //        for (size_t ii = 0; ii < dimRange; ++ii) {
      //          const CGJacobianType jacobian_ii(DomainType(-1),
      //                                           DomainType(1),
      //                                           DSC::FieldVector< size_t, dimDomain >(velocity_grid_view_->size(0)),
      //                                           basefunction_jacobians_values[ii]);
      //          basefunction_jacobians.emplace_back(jacobian_ii);
      //        }

      //        // calculate matrices
      //        VelocityExpressionFunctionType v("v", "v[0]", 1);
      //        VelocityExpressionFunctionType onebeam_left_boundary("v", "3*exp(3*v[0]+3)/(exp(6)-1)", 10);
      //        VelocityExpressionFunctionType one_minus_v_squared("v", "1-(v[0]^2)", 2);
      //        const typename Dune::GDT::Products::L2< VelocityGridViewType, R > l2_product(*velocity_grid_view_);
      //        for (size_t ii = 0; ii < dimRange; ++ii) {
      //          // Note: this assumes basefunctions_[0] is the constant function with value 1!!
      //          base_integrated_[ii] = l2_product.apply2(basefunctions_[0], basefunctions_[ii]);
      //          onebeam_left_boundary_values_[ii] = l2_product.apply2(onebeam_left_boundary, basefunctions_[ii]);
      //          for (size_t jj = 0; jj < dimRange; ++jj) {
      //            M_[ii][jj] = l2_product.apply2(basefunctions_[jj], basefunctions_[ii]);
      //            const auto v_times_base = DS::Functions::Product< VelocityFunctionType,
      //                                                              CGFunctionType >(v, basefunctions_[jj]);
      //            const auto jacobian_times_one_minus_v_squared
      //                = DS::Functions::Product< VelocityFunctionType, CGJacobianType >(one_minus_v_squared,
      //                                                                                 basefunction_jacobians[jj]);
      //            D_[ii][jj] = l2_product.apply2(v_times_base, basefunctions_[ii]);
      //            S_[ii][jj] = l2_product.apply2(jacobian_times_one_minus_v_squared, basefunction_jacobians[ii]);
      //          }
      //        }
      //        M_inverse_ = M_;
      //        M_inverse_.invert();
      //        is_calculated_ = true;
      //      }
    } // ... calculate()

    static MatrixType M_, D_, S_, M_inverse_;
    static RangeType base_integrated_;
    static bool is_calculated_;
    static RangeType basefunctions_values_at_minusone_;
    static RangeType basefunctions_values_at_plusone_;
    static std::string filename_;
    static bool exact_legendre_;
    static RangeType onebeam_left_boundary_values_;
    //    static std::vector< CGFunctionType > basefunctions_;
    //    static std::shared_ptr< const VelocityGridType > velocity_grid_;
    //    static std::shared_ptr< const VelocityGridViewType > velocity_grid_view_;
  };

public:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0]";
    grid_config["upper_right"]  = "[1.0]";
    grid_config["num_elements"] = "[100]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  static std::unique_ptr<ThisType> create(const ConfigType cfg = default_config(),
                                          const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr<const DefaultFluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const DefaultRHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config   = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  /** Reads basefunctions from a file instead of using the Legendre polynomials. Only CG functions with polOrder 1 are
   * supported. The *.csv file at the path basefunctions_file must contain the DoFs for the velocity grid in the
   * first row (e.g. -1.0, 0.0, 1.0 for a 2-cell grid with cells [-1,0] and [0,1]) and in each following row the
   * corresponding values of a basefunction at the DoFs ( e.g. 1, 1, 1 in the second row for the zero order Legendre
   * polynomial and (-1, 0, 1) in the third row for the first order Legendre polynomial ...).
   * */
  static std::unique_ptr<ThisType> create(const std::string basefunctions_file)
  {
    return create(default_config(basefunctions_file), static_id());
  }

  static ConfigType default_config(const std::string basefunctions_file = "", const std::string sub_name = "")
  {
    if (!(basefunctions_file.empty())) {
      GetData::set_filename(basefunctions_file);
      GetData::set_exact_legendre(false);
    }
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    flux_config["A"]    = GetData::create_flux_matrix();
    flux_config["b"] = DSC::to_string(RangeType(0));
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"]   = "[0.0]";
    rhs_config["upper_right"]  = "[1.0]";
    rhs_config["num_elements"] = "[1]";
    GetData::create_rhs_values(rhs_config);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"]   = "[0.0]";
    initial_value_config["upper_right"]  = "[1.0]";
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"]     = "x";
    initial_value_config["values.0"]     = GetData::create_initial_values();
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config;
    boundary_value_config["type"]       = BoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = GetData::create_boundary_values();
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  TwoBeams(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> rhs_in,
           const std::shared_ptr<const InitialValueType> initial_values_in, const ConfigType& grid_config_in,
           const ConfigType& boundary_info_in, const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, rhs_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }

  virtual double CFL() const override
  {
    return 0.4;
  }

  virtual double t_end() const override
  {
    return 4.0;
  }

  virtual bool is_linear() const override
  {
    return true;
  }

  virtual bool has_non_zero_rhs() const override
  {
    return true;
  }
}; // ... TwoBeams ...

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::MatrixType TwoBeams<E, D, d, R, rangeDim>::GetData::M_;

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::MatrixType TwoBeams<E, D, d, R, rangeDim>::GetData::D_;

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::MatrixType TwoBeams<E, D, d, R, rangeDim>::GetData::S_;

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::MatrixType TwoBeams<E, D, d, R, rangeDim>::GetData::M_inverse_;

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::RangeType TwoBeams<E, D, d, R, rangeDim>::GetData::base_integrated_(rangeDim);

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::RangeType
    TwoBeams<E, D, d, R, rangeDim>::GetData::basefunctions_values_at_minusone_(rangeDim);

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::RangeType
    TwoBeams<E, D, d, R, rangeDim>::GetData::basefunctions_values_at_plusone_(rangeDim);

template <class E, class D, size_t d, class R, size_t rangeDim>
bool TwoBeams<E, D, d, R, rangeDim>::GetData::is_calculated_(false);

template <class E, class D, size_t d, class R, size_t rangeDim>
std::string TwoBeams<E, D, d, R, rangeDim>::GetData::filename_("");

template <class E, class D, size_t d, class R, size_t rangeDim>
bool TwoBeams<E, D, d, R, rangeDim>::GetData::exact_legendre_(true);

template <class E, class D, size_t d, class R, size_t rangeDim>
typename TwoBeams<E, D, d, R, rangeDim>::RangeType
    TwoBeams<E, D, d, R, rangeDim>::GetData::onebeam_left_boundary_values_(rangeDim);

// template< class E, class D, size_t d, class R, size_t rangeDim >
// std::vector< typename TwoBeams< E, D, d, R, rangeDim >::GetData::CGFunctionType >
// TwoBeams< E, D, d, R, rangeDim >::GetData::basefunctions_;

// template< class E, class D, size_t d, class R, size_t rangeDim >
// std::shared_ptr< const typename TwoBeams< E, D, d, R, rangeDim >::GetData::VelocityGridViewType >
// TwoBeams< E, D, d, R, rangeDim >::GetData::velocity_grid_view_;

// template< class E, class D, size_t d, class R, size_t rangeDim >
// std::shared_ptr< const typename TwoBeams< E, D, d, R, rangeDim >::GetData::VelocityGridType >
// TwoBeams< E, D, d, R, rangeDim >::GetData::velocity_grid_;

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
