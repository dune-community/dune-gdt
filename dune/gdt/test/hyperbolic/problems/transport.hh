// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/stuff/functions/composition.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/test/nonstationary-eocstudy.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {

template< class EntityImp, class DomainFieldImp, size_t domainDim >
class PeriodicTransportFunction
  : public DS::GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, DomainFieldImp, domainDim, 1 >
{
  typedef DS::GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, DomainFieldImp, domainDim, 1 > BaseType;
  typedef PeriodicTransportFunction< EntityImp, DomainFieldImp, domainDim >      ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".periodictransport";
  }

  explicit PeriodicTransportFunction(const DomainType velocity,
                                     const double t,
                                     const DomainType lower_left,
                                     const DomainType upper_right,
                                     const std::string nm = static_id())
    : velocity_(velocity)
    , t_(t)
    , lower_left_(lower_left)
    , upper_right_(upper_right)
    , name_(nm)
  {}

  PeriodicTransportFunction(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".periodictransport";
  }

  virtual size_t order() const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    for (size_t ii = 0; ii < dimRange; ++ii) {
      ret[ii] = x[ii] - velocity_[ii]*t_;
      if (ret[ii] < lower_left_[ii] || ret[ii] > upper_right_[ii])
        ret[ii] = ret[ii] - std::floor(ret[ii]);
    }
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const override final
  {
    ret = JacobianRangeType(1);
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  const DomainType velocity_;
  const double t_;
  const DomainType lower_left_;
  const DomainType upper_right_;
  const std::string name_;
};

template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols >
class InitialValues
    : public DS::GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef typename DS::GlobalFunctionInterface< EntityImp,
                                                DomainFieldImp, domainDim,
                                                RangeFieldImp, rangeDim, rangeDimCols >     BaseType;
public:
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  InitialValues()
  {}

  virtual size_t order() const override
  {
   return 2;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const override
  {
    evaluate_helper(xx, ret, DS::Functions::internal::ChooseVariant< dimDomain >());
  }

private:
  void evaluate_helper(const DomainType& xx, RangeType& ret, const DS::Functions::internal::ChooseVariant< 1 >) const
  {
    if (DSC::FloatCmp::ge(xx[0], 0.2) && xx[0] < 0.4)
      ret[0] = 10000*std::pow(xx[0]-0.2,2)*std::pow(xx[0]-0.4,2)*std::exp(0.02-std::pow(xx[0]-0.2,2)-std::pow(xx[0]-0.4,2));
    else if (DSC::FloatCmp::ge(xx[0], 0.6) && xx[0] < 0.8)
      ret[0] = 1;
    else
      ret[0] = 0;
  }

  void evaluate_helper(const DomainType& xx, RangeType& ret, const DS::Functions::internal::ChooseVariant< 2 >) const
  {
    if (DSC::FloatCmp::ge(xx[0], 0.2) && xx[0] < 0.4 && DSC::FloatCmp::ge(xx[1], 0.2) && xx[1] < 0.4)
      ret[0] = 10000*std::pow(xx[0]-0.2,2)*std::pow(xx[0]-0.4,2)*std::exp(0.02-std::pow(xx[0]-0.2,2)-std::pow(xx[0]-0.4,2))*10000*std::pow(xx[1]-0.2,2)*std::pow(xx[1]-0.4,2)*std::exp(0.02-std::pow(xx[1]-0.2,2)-std::pow(xx[1]-0.4,2));
    else if (DSC::FloatCmp::ge(xx[0], 0.6) && xx[0] < 0.8 && DSC::FloatCmp::ge(xx[1], 0.6) && xx[1] < 0.8)
      ret[0] = 1;
    else
      ret[0] = 0;
  }
};


template< class LocalizableFunctionType, class GridViewType >
class TransportSolution
    : public DS::TimeDependentFunctionInterface
                      < typename DS::LocalizableFunctionInterface< typename LocalizableFunctionType::EntityType,
                                                                   typename LocalizableFunctionType::DomainFieldType,
                                                                   LocalizableFunctionType::dimDomain,
                                                                   typename LocalizableFunctionType::RangeFieldType,
                                                                   LocalizableFunctionType::dimRange,
                                                                   LocalizableFunctionType::dimRangeCols >,
                        double >
{
  typedef typename DS::TimeDependentFunctionInterface
                            < typename DS::LocalizableFunctionInterface< typename LocalizableFunctionType::EntityType,
                                                                         typename LocalizableFunctionType::DomainFieldType,
                                                                         LocalizableFunctionType::dimDomain,
                                                                         typename LocalizableFunctionType::RangeFieldType,
                                                                         LocalizableFunctionType::dimRange,
                                                                         LocalizableFunctionType::dimRangeCols >,
                              double >                                                                      BaseType;
  using typename BaseType::TimeIndependentFunctionType;
  typedef PeriodicTransportFunction< typename LocalizableFunctionType::EntityType,
                                     typename LocalizableFunctionType::DomainFieldType,
                                     LocalizableFunctionType::dimDomain >                   DomainTransportFunctionType;

  typedef typename DomainTransportFunctionType::DomainType DomainType;

public:
  TransportSolution(const LocalizableFunctionType localizable_func,
                    const GridViewType& grid_view,
                    const DomainType velocity,
                    const DomainType lower_left,
                    const DomainType upper_right)
    : localizable_func_(localizable_func)
    , grid_view_(grid_view)
    , velocity_(velocity)
    , lower_left_(lower_left)
    , upper_right_(upper_right)
  {}

  virtual std::unique_ptr< TimeIndependentFunctionType > evaluate_at_time(const double t) const
  {
    DomainTransportFunctionType x_minus_t(velocity_, t, lower_left_, upper_right_);
    typedef typename DS::Functions::Composition< DomainTransportFunctionType,
                                                 LocalizableFunctionType,
                                                 GridViewType >                 CompositionType;
    return DSC::make_unique< CompositionType >(x_minus_t, localizable_func_, grid_view_);
  }

  virtual std::string type() const
  {
    return "gdt.transportsolution";
  }

  virtual std::string name() const
  {
    return "gdt.transportsolution";
  }

private:
  LocalizableFunctionType localizable_func_;
  const GridViewType& grid_view_;
  const DomainType velocity_;
  const DomainType lower_left_;
  const DomainType upper_right_;
};


namespace Problems {


template< class E, class D, size_t d, class R, size_t r, size_t rC = 1 >
class Transport
  : public Default< E, D, d, R, r, rC >
{
  typedef Transport< E, D, d, R, r, rC > ThisType;
  typedef Default< E, D, d, R, r, rC > BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DummyEntityType;
  typedef typename Dune::Stuff::Functions::Affine< DummyEntityType,
                                                   R,
                                                   dimRange,
                                                   R,
                                                   dimRange,
                                                   dimDomain >                      FluxAffineFunctionType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux< FluxAffineFunctionType, E, D, d, R, r, rC > DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".transport";
  }

  virtual std::string type() const override
  {
    return BaseType::type() + ".transport";
  }

  static std::string short_id()
  {
    return "Transport";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0 0.0 0.0]";
    grid_config["upper_right"] = "[1.0 1.0 1.0]";
    grid_config["num_elements"] = "[8 8 8]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "periodic";
    return boundary_config;
  }

  static std::unique_ptr< ThisType > create(const ConfigType cfg = default_config(),
                                            const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr< const DefaultFluxType > flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr< const DefaultRHSType > rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr< const DefaultInitialValueType > initial_values(DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr< const DefaultBoundaryValueType > boundary_values(DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique< ThisType >(flux, rhs, initial_values,
                                                  grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config();
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType flux_config;
    flux_config["A.0"] = "[1]";
    flux_config["A.1"] = "[2]";
    flux_config["b"] = "[0 0; 0 0]";
    config.add(flux_config, "flux", true);
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0 0.0 0.0]";
    initial_value_config["upper_right"] = "[1.0 1.0 1.0]";
    if (dimDomain == 1)
      initial_value_config["num_elements"] = "[5]";
    else
      initial_value_config["num_elements"] = "[5 5 1]";
    initial_value_config["variable"] = "x";
    if (dimDomain == 1)
      initial_value_config["values"] = "[0.0 10000*((x[0]-0.2)^2)*((x[0]-0.4)^2)*exp(0.02-((x[0]-0.2)^2)-((x[0]-0.4)^2)) 0.0 1.0 0.0]"; //"[0 sin(pi/2+5*pi*(x[0]-0.3))*exp(-(200*(x[0]-0.3)*(x[0]-0.3))) 0 1.0 0.0]";
    else
      initial_value_config["values"] = std::string("[0 0 0 0 0 ") +
                                       std::string( "0 10000*((x[0]-0.2)^2)*((x[0]-0.4)^2)*exp(0.02-((x[0]-0.2)^2)-((x[0]-0.4)^2))*10000*((x[1]-0.2)^2)*((x[1]-0.4)^2)*exp(0.02-((x[1]-0.2)^2)-((x[1]-0.4)^2)) 0 0 0 ") +
                                       std::string( "0 0 0 0 0 ") +
                                       std::string( "0 0 0 1 0 ") +
                                       std::string( "0 0 0 0 0]");
    initial_value_config["order"] = "10";
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values", true);
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  Transport(const std::shared_ptr< const FluxType > flux_in,
            const std::shared_ptr< const RHSType > rhs_in,
            const std::shared_ptr< const InitialValueType > initial_values_in,
            const ConfigType& grid_config_in,
            const ConfigType& boundary_info_in,
            const std::shared_ptr< const BoundaryValueType > boundary_values_in)
    : BaseType(flux_in,
               rhs_in,
               initial_values_in,
               grid_config_in,
               boundary_info_in,
               boundary_values_in)
  {}

  virtual double CFL() const override
  {
    if (dimDomain == 1)
      return 0.5;
    else
      return 0.3;
  }

  virtual double t_end() const override
  {
    return 1.0;
  }

  virtual bool is_linear() const override
  {
    return true;
  }
};

} // namespace Problems

template< class G, class R = double, size_t r = 1, size_t rC = 1 >
class TransportTestCase
  : public Dune::GDT::Tests::NonStationaryTestCase< G, Problems::Transport< typename G::template Codim< 0 >::Entity,
                                                                     typename G::ctype, G::dimension, R, r, rC > >
{
  typedef typename G::template Codim< 0 >::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
public:
  static const size_t dimRange = r;
  static const size_t dimRangeCols = rC;
  typedef typename Problems::Transport< E, D, d, R, r, rC > ProblemType;
private:
  typedef typename Dune::GDT::Tests::NonStationaryTestCase< G, ProblemType > BaseType;
public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  TransportTestCase(const size_t num_refs = (d == 1 ? 4 : 2))
    : BaseType(Stuff::Grid::Providers::Cube< G >::create(ProblemType::default_grid_config())->grid_ptr(), num_refs)
    , reference_grid_view_(BaseType::reference_grid_view())
    , problem_(*(ProblemType::create(ProblemType::default_config())))
  {
    typedef InitialValues< E, D, d, R, r, 1 > LocalizableInitialValueType;
    const LocalizableInitialValueType initial_values;
    exact_solution_ = std::make_shared< const TransportSolution< LocalizableInitialValueType,
                                                           LevelGridViewType > >(initial_values,
                                                                                 reference_grid_view_,
                                                                                 DSC::from_string< typename DSC::FieldVector< D, d > >("[1.0 2.0]"),
                                                                                 DSC::from_string< typename DSC::FieldVector< D, d > >(problem_.grid_config()["lower_left"]),
                                                                                 DSC::from_string< typename DSC::FieldVector< D, d > >(problem_.grid_config()["upper_right"]));
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return true;
  }

  virtual std::bitset< d > periodic_directions() const override final
  {
    std::bitset< d > periodic_dirs;
    periodic_dirs.set();
    return periodic_dirs;
  }

  virtual const std::shared_ptr< const SolutionType > exact_solution() const override final
  {
    return std::shared_ptr< const SolutionType >(exact_solution_);
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    std::string domainstring;
    if (d == 1)
      domainstring = "||  domain = [0, 1]                                                   ||\n";
    else
      domainstring = "||  domain = [0, 1] x [0, 1]                                          ||\n";
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: Transport                                               ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        <<    domainstring
        << "||  flux = u[0]                                                       ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: exact solution                                ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const LevelGridViewType reference_grid_view_;
  const ProblemType problem_;
  std::shared_ptr< const SolutionType > exact_solution_;
}; // class TransportTestCase

} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH
