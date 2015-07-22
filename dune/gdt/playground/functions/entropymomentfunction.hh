// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_FUNCTIONS_ENTROPYMOMENT_HH
#define DUNE_STUFF_FUNCTIONS_ENTROPYMOMENT_HH

#include <memory>
#include <cmath>
#include <algorithm>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/functions/constant.hh>

#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/products/l2.hh>

#include "dune/stuff/functions/interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


template< class LocalizableFunctionTypeLeft, class LocalizableFunctionTypeRight >
class DotProduct
  : public LocalizableFunctionInterface< typename LocalizableFunctionTypeLeft::EntityType,
                                         typename LocalizableFunctionTypeLeft::DomainFieldType,
                                         LocalizableFunctionTypeLeft::dimDomain,
                                         typename LocalizableFunctionTypeLeft::RangeFieldType,
                                         1,
                                         1 >
{
  typedef LocalizableFunctionInterface< typename LocalizableFunctionTypeLeft::EntityType,
                                        typename LocalizableFunctionTypeLeft::DomainFieldType,
                                        LocalizableFunctionTypeLeft::dimDomain,
                                        typename LocalizableFunctionTypeLeft::RangeFieldType,
                                        1,
                                        1 >                                                       BaseType;
  typedef DotProduct< LocalizableFunctionTypeLeft, LocalizableFunctionTypeRight >                 ThisType;
public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;

private:
  class Localfunction
    : public LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >
  {
    typedef LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols > BaseType;
  public:
    Localfunction(const EntityType& entity,
                  const std::vector< LocalizableFunctionTypeLeft >& localizable_function_one,
                  const std::vector< LocalizableFunctionTypeRight >& localizable_function_two)
      : BaseType(entity)
      , entity_(entity)
      , localizable_function_one_(localizable_function_one)
      , localizable_function_two_(localizable_function_two)
    {}

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return localizable_function_one_[0].local_function(entity_)->order()
           + localizable_function_two_[0].local_function(entity_)->order();
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      ret = RangeType(0);
      for (size_t ii = 0; ii < localizable_function_one_.size(); ++ii) {
        ret += localizable_function_one_[ii].local_function(entity_)->evaluate(xx) * localizable_function_two_[ii].local_function(entity_)->evaluate(xx);
      }
    }

    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const override
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const EntityType& entity_;
    const std::vector< LocalizableFunctionTypeLeft >& localizable_function_one_;
    const std::vector< LocalizableFunctionTypeRight >& localizable_function_two_;
  }; // class Localfunction

public:
  using typename BaseType::LocalfunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".composition";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["values"] = "[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
//    // get correct config
//    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
//    const Common::Configuration default_cfg = default_config();
//    // calculate number of values and get values
//    auto num_elements = cfg.get("num_elements",
//                                default_cfg.get< Common::FieldVector< size_t, dimDomain > >("num_elements"), dimDomain);
//    size_t num_values = 1;
//    for (size_t ii = 0; ii < num_elements.size(); ++ii)
//      num_values *= num_elements[ii];
//    std::vector< RangeType > values(num_values);
//    if (config.has_key("values.0")) { // get every value from its own config entry
//      try { // get value directly as RangeType
//        for (size_t ii = 0; ii < num_values; ++ii)
//          values[ii] = cfg.get< RangeType >("values." + DSC::toString(ii),
//                                            dimRange,
//                                            dimRangeCols);
//      } catch (const Exceptions::conversion_error& e) {
//        if (dimRangeCols == 1) { // get every value from its own config entry as the first col of a matrix
//          for (size_t ii = 0; ii < num_values; ++ii) {
//            const auto values_matrix = cfg.get< Common::FieldMatrix< RangeFieldType, dimRange, 1 > >("values." + DSC::toString(ii),
//                                                                                                     dimRange,
//                                                                                                     1);
//            // this fromString(toString(...)) construct avoids a compilation error if dimRangeCols > 1 and is easier
//            // than creating templated helper methods
//            values[ii] = DSC::fromString< RangeType >(DSC::toString(values_matrix[ii]));
//          }
//        } else {
//          std::cout << e.what() << std::endl;
//        }
//      }
//    } else {
//      // get values as a vector of scalars
//      auto values_rf = cfg.get("values", default_cfg.get< std::vector< RangeFieldType > >("values"), num_values);
//      for (size_t ii = 0; ii < values_rf.size(); ++ii)
//        values[ii] = RangeType(values_rf[ii]);
//    }
//    // create
//    return Common::make_unique< ThisType >(
//            cfg.get("lower_left",
//                    default_cfg.get< Common::FieldVector< DomainFieldType, dimDomain > >("lower_left"), dimDomain),
//            cfg.get("upper_right",
//                    default_cfg.get< Common::FieldVector< DomainFieldType, dimDomain > >("upper_right"), dimDomain),
//            std::move(num_elements),
//            std::move(values),
//            cfg.get("name", default_cfg.get< std::string > ("name")));
  } // ... create(...)

  DotProduct(const std::vector< LocalizableFunctionTypeLeft > localizable_functions_one,
             const std::vector< LocalizableFunctionTypeRight > localizable_functions_two,
             const std::string nm = static_id())
    : localizable_functions_one_(localizable_functions_one)
    , localizable_functions_two_(localizable_functions_two)
    , name_(nm)
  {}

  DotProduct(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".composition";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& entity) const override
  {
     return DSC::make_unique< Localfunction >(entity,
                                              localizable_functions_one_,
                                              localizable_functions_two_);
  } // ... local_function(...)

private:
  std::vector< LocalizableFunctionTypeLeft > localizable_functions_one_;
  std::vector< LocalizableFunctionTypeRight > localizable_functions_two_;
  std::string name_;
}; // class DotProduct


/**
 * \brief Simple affine function of the form f(x) = A*x + b
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class EntropyMomentFlux
{
  EntropyMomentFlux() { static_assert(AlwaysFalse< EntityImp >::value, "Not available for rangeDimCols > 1!"); }
};


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class EntropyMomentFlux< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > BaseType;
  typedef EntropyMomentFlux< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >       ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  // these should be template arguments?
  typedef double                                                                    VelocityFieldImp;
  typedef typename Dune::SGrid< 1, 1, VelocityFieldImp >            VelocityGridType;
  typedef Dune::Stuff::Grid::Providers::Cube< VelocityGridType >                    VelocityGridProviderType;
  typedef typename VelocityGridType::LeafGridView                                   VelocityGridViewType;
  typedef typename VelocityGridType::template Codim< 0 >::Entity                    VelocityEntityType;
  typedef typename DS::LocalizableFunctionInterface< VelocityEntityType,
                                                     VelocityFieldImp, 1,
                                                     RangeFieldImp, 1, 1 >          VelocityFunctionType;
  typedef typename DS::Functions::Expression< VelocityEntityType,
                                              VelocityFieldImp, 1,
                                              RangeFieldImp, 1, 1 >                 VelocityExpressionFunctionType;

  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldImp >              VectorType;
  typedef typename Dune::GDT::Spaces::CGProvider< VelocityGridType,
                                                  DSG::ChooseLayer::leaf,
                                                  Dune::GDT::ChooseSpaceBackend::pdelab,
                                                  1, RangeFieldImp, 1, 1 >          CGProviderType;
  typedef typename CGProviderType::Type                                             CGSpaceType;
  typedef Dune::GDT::DiscreteFunction< CGSpaceType, VectorType >                    CGFunctionType;
  typedef typename DS::Functions::Constant< VelocityEntityType, VelocityFieldImp, 1, RangeFieldImp, 1 > ConstantFunctionType;
  typedef typename Dune::FieldMatrix< RangeFieldImp, rangeDim, domainDim >  MatrixType;

  using typename BaseType::LocalfunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".affine";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["A"] = internal::Get< RangeFieldImp, rangeDim, domainDim >::value_str();
    config["b"] = internal::Get< RangeFieldImp, rangeDim, 1 >::value_str();
    config["sparse"] = "false";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
//    // get correct config
//    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
//    const Common::Configuration default_cfg = default_config();
//    return Common::make_unique< ThisType >(
//          cfg.get("A", default_cfg.get< MatrixType >("A")),
//          cfg.get("b", default_cfg.get< RangeType >("b")),
//          cfg.get("sparse", default_cfg.get< bool >("sparse")),
//          cfg.get("name",  default_cfg.get< std::string >("name")));
  } // ... create(...)

  explicit EntropyMomentFlux(std::shared_ptr< const VelocityGridViewType > velocity_grid_view,
                             const std::vector< CGFunctionType >& basefunctions,
                             const RangeType alpha,
                             const RangeFieldImp beta,
                             const RangeFieldImp tau,
                             const RangeFieldImp epsilon_gamma,
                             const RangeFieldImp xi,
                             const std::string name = static_id())
    : velocity_grid_view_(velocity_grid_view)
    , basefunctions_(basefunctions)
    , m_m_T_(rangeDim)
    , name_(name)
    , alpha_(alpha)
    , v_("v", "v[0]", 1)
    , exp_("v", "exp(v[0])", 40)
    , v_m_()
    , beta_(beta)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , zeta_(0)
    , xi_(xi)
    , l2_product_(*velocity_grid_view_)
  {
    const auto it_end = velocity_grid_view_->template end< 0 >();
    for (auto it = velocity_grid_view_->template begin< 0 >(); it != it_end; ++it) {
      auto& entity = *it;
      const auto intersection_it_end = velocity_grid_view_->iend(entity);
      for (auto intersection_it = velocity_grid_view_->ibegin(entity); intersection_it != intersection_it_end; ++intersection_it)
      {
        auto& intersection = *intersection_it;
        const auto intersection_center = intersection.geometryInInside().center();
        RangeType base_evaluated;
        for (size_t ii = 0; ii < rangeDim; ++ii) {
          base_evaluated[ii] = basefunctions_[ii].local_function(entity)->evaluate(intersection_center);
        }
        const auto norm_m = base_evaluated.two_norm();
        if (norm_m > zeta_)
          zeta_ = norm_m;
      }
    }
    std::cout << "zeta: " << zeta_  << std::endl;
    for (size_t ii = 0; ii < rangeDim; ++ii) {
      v_m_.emplace_back(typename DS::Functions::Product< VelocityFunctionType, CGFunctionType >(v_, basefunctions_[ii]));
      for (size_t jj = 0; jj < domainDim; ++jj) {
         m_m_T_[ii].emplace_back(typename DS::Functions::Product< CGFunctionType, CGFunctionType >(basefunctions_[ii],
                                                                                                   basefunctions_[jj]));
      }
    }
        alpha_ = RangeType(0);
        alpha_[0] = std::log(2);
  }

  EntropyMomentFlux(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".affine";
  }

  virtual size_t order() const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    for (size_t ii = 0; ii < last_four_x_.size(); ++ii) {
      if (last_four_x_[ii] == x) {
        ret = last_four_alpha_[ii];
        return;
      }
    }
    current_x_ = x;
    last_four_x_.emplace(last_four_x_.begin(), x);
    if (last_four_x_.size() > 4)
      last_four_x_.pop_back();
    evaluate(ret);
  }

  void evaluate(RangeType& ret) const
  {
    // reset alpha_
    alpha_ = RangeType(0);
    alpha_[0] = std::log(2);
    size_t quadrature_order = 1;
    // reset r_
    r_ = 1e-8;
    ignore_order_ = false;

    // compute f_alpha
    RangeFieldImp current_f_alpha;
    RangeType current_g_alpha, current_d_alpha;
    compute_data(alpha_, quadrature_order, current_f_alpha, current_g_alpha, current_d_alpha);

    auto g_alpha_norm = current_g_alpha.two_norm();
    auto d_alpha_norm = current_d_alpha.two_norm();
    auto condition_2 = 5.0*zeta_*d_alpha_norm;

//    std::cout << "initial conditions: " << std::endl;
//    std::cout << " alpha: " << DSC::toString(alpha_);
//    std::cout << " x: " << DSC::toString(current_x);
//    std::cout << " H_alpha: " << DSC::toString(H_alpha) << std::endl;
//    std::cout << "g_alpha_ " << DSC::toString(g_alpha);
//    std::cout << " d_alpha " << DSC::toString(d_alpha) << std::endl;

    std::vector< RangeType > alpha_iterates(1, alpha_);
    while (g_alpha_norm > tau_ && condition_2 > std::log(1.0 + epsilon_gamma_)) {
      // find zero of g_alpha by Newton Method with Armijo backtracking

      // backtracking line search
      RangeFieldImp factor = 1;
      while (quadrature_order <= 50) {
        auto new_alpha = armijo_backtracking(alpha_, quadrature_order, current_f_alpha, current_g_alpha, current_d_alpha, factor);
        // compute data with higher quadrature order ...
        RangeFieldImp current_f_alpha_finer;
        RangeType current_g_alpha_finer, current_d_alpha_finer;
        compute_data(alpha_, quadrature_order + 10, current_f_alpha_finer, current_g_alpha_finer, current_d_alpha_finer);
        auto new_alpha_finer = current_d_alpha_finer;
        new_alpha_finer *= factor;
        new_alpha_finer += alpha_;
        auto new_f_alpha_finer = compute_f_alpha(new_alpha_finer, quadrature_order + 10);
        // ... and check if backtracking condition is still fulfilled
        if (std::abs(current_f_alpha_finer - current_f_alpha) > 0.5e-12 || new_f_alpha_finer > current_f_alpha_finer + current_g_alpha_finer*current_d_alpha_finer*xi_*factor) {
          // go on with the backtracking using the finer quadrature
          alpha_ = new_alpha_finer;
          quadrature_order += 1;
          if (quadrature_order > 11 && !ignore_order_) {
            alpha_ = RangeType(0);
            alpha_[0] = std::log(2);
            realizability_projection();
            quadrature_order = 1;
            std::cout << "Das hätte ich ja nicht gedacht" << std::endl;
            compute_data(alpha_, quadrature_order, current_f_alpha, current_g_alpha, current_d_alpha);
            break;
          } else {
            compute_data(alpha_, quadrature_order, current_f_alpha, current_g_alpha, current_d_alpha);
          }
        } else {
          alpha_ = new_alpha;
          alpha_iterates.emplace_back(alpha_);
          compute_data(alpha_, quadrature_order, current_f_alpha, current_g_alpha, current_d_alpha);
          break;
        }
      }

      // compute loop conditions
      g_alpha_norm = current_g_alpha.two_norm();
      d_alpha_norm = current_d_alpha.two_norm();
      condition_2 = 5.0*zeta_*d_alpha_norm;
//            std::cout << DSC::toString(g_alpha_norm) << " AND " << DSC::toString(condition_2) << std::endl;
      //      std::cout << "alpha: " << DSC::toString(alpha_) << std::endl;
    }

    // calculate f(u) = < v m G_alpha_hat >
    for (size_t ii = 0; ii < rangeDim; ++ii) {
      ret[ii] = l2_product_.apply2(v_m_[ii], *G_alpha_);
    }
    last_four_alpha_.emplace(last_four_alpha_.begin(), ret);
    if (last_four_alpha_.size() > 4)
      last_four_alpha_.pop_back();
    std::cout << "Evaluation done" << std::endl;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const override final
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  void compute_data(RangeType& alpha, const size_t order, RangeFieldImp& f_alpha, RangeType& g_alpha, RangeType& d_alpha, const bool calculate_f_alpha = true) const
  {
    if (calculate_f_alpha)
      f_alpha = compute_f_alpha(alpha, order);

    // calculate g(\alpha) = < m G_\alpha > - x ...
    g_alpha = current_x_;
    g_alpha *= -1.0;
    // ... and H_(\alpha) = < m m^T G_\alpha > ...
    // ... by using the L2 Product on the velocity grid view
    for (size_t ii = 0; ii < rangeDim; ++ii) {
      const auto base_ii_times_G_alpha = DS::Functions::Product< CGFunctionType, Composition< VelocityExpressionFunctionType,
          DotProduct< ConstantFunctionType, CGFunctionType > > >(basefunctions_[ii], *G_alpha_);
      g_alpha[ii] += quadrature(base_ii_times_G_alpha, order);
      for (size_t jj = 0; jj < rangeDim; ++jj) {
        const auto m_m_T_ii_jj_times_G_alpha
            = DS::Functions::Product< typename DS::Functions::Product< CGFunctionType, CGFunctionType >,
                                      Composition< VelocityExpressionFunctionType,
                                                   DotProduct< ConstantFunctionType, CGFunctionType > > >(m_m_T_[ii][jj], *G_alpha_);
        H_alpha_[ii][jj] = quadrature(m_m_T_ii_jj_times_G_alpha, order);
      }
    }
      H_alpha_.invert();
      H_alpha_ *= -1.0;
      H_alpha_.mv(g_alpha, d_alpha);



//    if (DSC::isnan(cond_H_alpha)
//        || DSC::isinf(cond_H_alpha)
//        || cond_H_alpha > 1e16)
//    {
//      alpha = RangeType(0);
//      alpha[0] = std::log(2);
//      realizability_projection();
//      compute_data(alpha, order, f_alpha, g_alpha, d_alpha);
//    }
  }

  RangeFieldImp compute_f_alpha(const RangeType& alpha, const size_t order) const
  {
    std::vector< ConstantFunctionType > alpha_function;
    for (size_t ii = 0; ii < rangeDim; ++ii) {
      alpha_function.emplace_back(ConstantFunctionType(alpha[ii]));
    }
    DotProduct< ConstantFunctionType, CGFunctionType > alpha_T_m = DotProduct< ConstantFunctionType, CGFunctionType >(alpha_function, basefunctions_);
    // f(\alpha) = < G_\alpha > - \alpha^T x
    G_alpha_ = std::make_shared< Composition< VelocityExpressionFunctionType,
                       DotProduct< ConstantFunctionType, CGFunctionType > > >(exp_, alpha_T_m);
    auto f_alpha = quadrature(*G_alpha_, order);
    for (size_t ii = 0; ii < rangeDim; ++ii)
       f_alpha -= alpha[ii]*current_x_[ii];
    return f_alpha;
  }

  void realizability_projection() const
  {
    // realizability projection x_new = (1-r)x + rQu
    DomainType x_new = current_x_;
    x_new *= 1 - r_;
    current_x_[0] *= r_;
    x_new[0] += current_x_[0];
    current_x_ = x_new;

    // increase value of r
    r_ += r_;
    if (r_ > 0.1) {
      ignore_order_ = true;
      std::cout << "r ist zu groß!!: " << DSC::toString(r_) << std::endl;
    }
  }

  RangeFieldImp quadrature(const VelocityFunctionType& integrand, const size_t order) const
  {
    RangeFieldImp ret(0);
    const auto it_end = velocity_grid_view_->template end< 0 >();
    for (auto it = velocity_grid_view_->template begin< 0 >(); it != it_end; ++it)
    {
      const auto& entity = *it;
      // quadrature
      const auto& volumeQuadrature = QuadratureRules< VelocityFieldImp, 1 >::rule(entity.type(),
                                                                                        boost::numeric_cast< int >(order));
      // loop over all quadrature points
      const auto local_function = integrand.local_function(entity);
      for (const auto& quadPoint : volumeQuadrature) {
        const auto x = quadPoint.position();
        // integration factors
        const auto integrationFactor = entity.geometry().integrationElement(x);
        const auto quadratureWeight = quadPoint.weight();
        // add to integral
        ret += local_function->evaluate(x)*integrationFactor*quadratureWeight;
      } // compute integral
    }
    return ret;
  }

  RangeType armijo_backtracking(const RangeType& alpha, const size_t order, const RangeFieldImp& f_alpha, const RangeType& g_alpha, const RangeType& direction, RangeFieldImp& factor) const
  {
    // calculate xi * g_alpha^T * step_length
    auto g_alpha_xi = g_alpha;
    g_alpha_xi *= xi_;
    auto step = direction;
    step *= factor;
    RangeType alpha_new = alpha + step;
    auto f_alpha_new = compute_f_alpha(alpha_new, order);
    auto armijo_condition = f_alpha + g_alpha_xi*direction;
    while (f_alpha_new > armijo_condition) {
      factor *= beta_;
      step = direction;
      step *= factor;
      alpha_new = alpha + step;
      f_alpha_new = compute_f_alpha(alpha_new, order);
      armijo_condition = f_alpha + g_alpha_xi*step;
    }
    return alpha_new;
  }

  std::shared_ptr< const VelocityGridViewType > velocity_grid_view_;
  const std::vector< CGFunctionType >& basefunctions_;
  std::vector< std::vector< typename DS::Functions::Product< CGFunctionType, CGFunctionType > > > m_m_T_;
  const std::string name_;
  mutable RangeType alpha_;
  const RangeFieldImp beta_;
  const RangeFieldImp tau_;
  const RangeFieldImp epsilon_gamma_;
  const RangeFieldImp xi_;
  const typename Dune::GDT::Products::L2< VelocityGridViewType, RangeFieldImp > l2_product_;
  mutable DomainType current_x_;
  mutable std::shared_ptr< Composition< VelocityExpressionFunctionType,
               DotProduct< ConstantFunctionType, CGFunctionType > > > G_alpha_;
  mutable JacobianRangeType H_alpha_;
  mutable RangeType d_alpha_;
  mutable double r_;
  RangeFieldImp zeta_;
  VelocityExpressionFunctionType v_;
  VelocityExpressionFunctionType exp_;
  std::vector< DS::Functions::Product< VelocityFunctionType, CGFunctionType > > v_m_;
  mutable std::vector< DomainType > last_four_x_;
  mutable std::vector< RangeType > last_four_alpha_;
  mutable bool ignore_order_;
};


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_ENTROPYMOMENT_HH
