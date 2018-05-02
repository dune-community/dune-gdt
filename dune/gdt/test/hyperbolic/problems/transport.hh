// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/composition.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace {


template <class DomainType>
double pow1(const DomainType& x, const size_t ii)
{
  return std::pow(x[ii] - 0.2, 2);
}
template <class DomainType>
double pow2(const DomainType& x, const size_t ii)
{
  return std::pow(x[ii] - 0.4, 2);
}
template <class DomainType>
double exp1(const DomainType& x, const size_t ii)
{
  return std::exp(0.02 - pow1(x, ii) - pow2(x, ii));
}

template <size_t dim>
struct initial_vals_helper
{
  template <class DomainType, class RangeType>
  static void evaluate(const DomainType& x, RangeType& ret)
  {
    if (Dune::XT::Common::FloatCmp::ge(x, DomainType(0.2)) && XT::Common::FloatCmp::lt(x, DomainType(0.4)))
      ret[0] = 10000 * pow1(x, 0) * pow2(x, 0) * exp1(x, 0);
    else if (Dune::XT::Common::FloatCmp::ge(x, DomainType(0.6)) && XT::Common::FloatCmp::lt(x, DomainType(0.8)))
      ret[0] = 1;
    else
      ret[0] = 0;
  }
}; // struct initial_vals_helper<1>

template <>
struct initial_vals_helper<2>
{
  template <class DomainType, class RangeType>
  static void evaluate(const DomainType& x, RangeType& ret)
  {
    if (Dune::XT::Common::FloatCmp::ge(x, DomainType(0.2)) && XT::Common::FloatCmp::lt(x, DomainType(0.4)))
      ret[0] = 10000 * pow1(x, 0) * pow2(x, 0) * exp1(x, 0) * 10000 * pow1(x, 1) * pow2(x, 1) * exp1(x, 1);
    else if (Dune::XT::Common::FloatCmp::ge(x, DomainType(0.6)) && XT::Common::FloatCmp::lt(x, DomainType(0.8)))
      ret[0] = 1;
    else
      ret[0] = 0;
  }
}; // struct initial_vals_helper<2>

template <>
struct initial_vals_helper<3>
{
  template <class DomainType, class RangeType>
  static void evaluate(const DomainType& x, RangeType& ret)
  {
    if (Dune::XT::Common::FloatCmp::ge(x, DomainType(0.2)) && XT::Common::FloatCmp::lt(x, DomainType(0.4)))
      ret[0] = 10000 * pow1(x, 0) * pow2(x, 0) * exp1(x, 0) * 10000 * pow1(x, 1) * pow2(x, 1) * exp1(x, 1) * 10000
               * pow1(x, 2) * pow2(x, 2) * exp1(x, 2);
    else if (Dune::XT::Common::FloatCmp::ge(x, DomainType(0.6)) && XT::Common::FloatCmp::lt(x, DomainType(0.8)))
      ret[0] = 1;
    else
      ret[0] = 0;
  }
}; // struct initial_vals_helper<3>


} // anonymous namespace


// A simple function x -> x - v t, where v is a velocity and t the current time.
// The range is restricted to a (multidimensional) interval [lower_left, upper_right] by applying periodic boundary
// conditions.
template <class EntityImp, class DomainFieldImp, size_t domainDim>
class PeriodicTransportFunction
    : public XT::Functions::GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, DomainFieldImp, domainDim, 1>
{
  typedef XT::Functions::GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, DomainFieldImp, domainDim, 1>
      BaseType;
  typedef PeriodicTransportFunction<EntityImp, DomainFieldImp, domainDim> ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
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
                                     const DomainType lower_left,
                                     const DomainType upper_right,
                                     const std::string nm = static_id())
    : velocity_(velocity)
    , lower_left_(lower_left)
    , upper_right_(upper_right)
    , name_(nm)
  {
  }

  PeriodicTransportFunction(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".periodictransport";
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret, const XT::Common::Parameter& mu = {}) const override final
  {

    const RangeFieldType t = mu.get("t")[0];
    for (size_t ii = 0; ii < dimRange; ++ii) {
      ret[ii] = x[ii] - velocity_[ii] * t;
      while (ret[ii] < lower_left_[ii])
        ret[ii] += upper_right_[ii] - lower_left_[ii];
      while (ret[ii] > upper_right_[ii])
        ret[ii] -= upper_right_[ii] - lower_left_[ii];
    }
  }

  virtual void jacobian(const DomainType& /*x*/,
                        JacobianRangeType& ret,
                        const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    ret = JacobianRangeType(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii][ii] = 1.;
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  const DomainType velocity_;
  const DomainType lower_left_;
  const DomainType upper_right_;
  const std::string name_;
};

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class TransportInitialValues
    : public XT::Functions::
          GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef typename XT::Functions::
      GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
          BaseType;

public:
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
  {
    return 50;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret, const XT::Common::Parameter& /*mu*/ = {}) const override
  {
    initial_vals_helper<dimDomain>::evaluate(x, ret);
  }
}; // class TransportInitialValues

template <class LocalizableFunctionType, class GridLayerType>
class TransportSolution
    : public XT::Functions::
          CompositionFunction<PeriodicTransportFunction<typename LocalizableFunctionType::EntityType,
                                                        typename LocalizableFunctionType::DomainFieldType,
                                                        LocalizableFunctionType::dimDomain>,
                              LocalizableFunctionType,
                              GridLayerType>
{
  typedef PeriodicTransportFunction<typename LocalizableFunctionType::EntityType,
                                    typename LocalizableFunctionType::DomainFieldType,
                                    LocalizableFunctionType::dimDomain>
      PeriodicTransportFunctionType;
  typedef typename PeriodicTransportFunctionType::DomainType DomainType;
  typedef XT::Functions::CompositionFunction<PeriodicTransportFunctionType, LocalizableFunctionType, GridLayerType>
      BaseType;

public:
  TransportSolution(const LocalizableFunctionType initial_values, const DomainType velocity)
    : BaseType(PeriodicTransportFunctionType(velocity, DomainType(0), DomainType(1)), initial_values)
  {
  }
}; // class TransportSolution<...>


namespace Problems {


template <class E, class D, size_t d, class U, class R, size_t r>
class Transport : public ProblemBase<E, D, d, U, R, r>
{
  typedef Transport<E, D, d, U, R, r> ThisType;
  typedef ProblemBase<E, D, d, U, R, r> BaseType;

public:
  static const bool linear = true;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  typedef typename XT::Functions::AffineFluxFunction<E, D, d, U, R, r, d> ActualFluxType;
  typedef typename XT::Functions::AffineFluxFunction<E, D, d, U, R, r, 1> ActualRhsType;
  typedef XT::Functions::GlobalLambdaFunction<E, D, d, R, r, 1> ActualBoundaryValueType;
  typedef ActualBoundaryValueType ActualInitialValueType;

  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;

  Transport(const XT::Common::Configuration& grid_cfg = BaseType::default_grid_cfg(),
            const XT::Common::Configuration& boundary_cfg = BaseType::default_boundary_cfg())
    : BaseType(create_flux(),
               create_rhs(),
               create_initial_values(grid_cfg),
               create_boundary_values(),
               grid_cfg,
               boundary_cfg,
               dimDomain == 1 ? 0.5 : (dimDomain == 2 ? 0.3 : 0.15),
               1.,
               false)
  {
  }

  static std::string static_id()
  {
    return "Transport";
  }

  static FluxType* create_flux()
  {
    FieldVector<MatrixType, dimDomain> A(MatrixType(0.));
    for (size_t ii = 0; ii < dimRange; ++ii) {
      A[0][ii][ii] = 1.;
      if (dimDomain > 1)
        A[1][ii][ii] = 2.;
      if (dimDomain > 2)
        A[2][ii][ii] = 3.;
    }
    return new ActualFluxType(A, typename ActualFluxType::RangeType(0));
  }

  static RhsType* create_rhs()
  {
    return new ActualRhsType(FieldVector<MatrixType, 1>(MatrixType(0.)));
  } // ... create_rhs(...)

  static InitialValueType* create_initial_values(const XT::Common::Configuration& /*grid_cfg*/)
  {
    return new ActualInitialValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          RangeType ret;
          initial_vals_helper<dimDomain>::evaluate(x, ret);
          return ret;
        },
        50);
  } // ... create_initial_values()

  virtual BoundaryValueType* create_boundary_values()
  {
    return new ActualBoundaryValueType([=](const DomainType&, const XT::Common::Parameter&) { return 0; }, 0);
  } // ... create_boundary_values()
}; // class Transport<...>


} // namespace Problems


template <class G, class R = double, size_t r = 1>
class TransportTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               Problems::Transport<
                                   typename G::template Codim<0>::Entity,
                                   typename G::ctype,
                                   G::dimension,
                                   typename internal::DiscreteFunctionProvider<G,
                                                                               GDT::SpaceType::product_fv,
                                                                               0,
                                                                               R,
                                                                               r,
                                                                               1,
                                                                               GDT::Backends::gdt>::type,
                                   R,
                                   r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  static const size_t dimRange = r;
  static const size_t dimRangeCols = 1;
  typedef
      typename internal::DiscreteFunctionProvider<G, GDT::SpaceType::product_fv, 0, R, r, 1, GDT::Backends::gdt>::type
          U;
  typedef typename Problems::Transport<E, D, d, U, R, r> ProblemType;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  TransportTestCase(const size_t num_refs = (d == 1 ? 4 : 2), const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
  {
    typedef TransportInitialValues<E, D, d, R, r, 1> LocalizableInitialValueType;
    exact_solution_ = std::make_shared<const TransportSolution<LocalizableInitialValueType, LevelGridViewType>>(
        LocalizableInitialValueType{},
        Dune::XT::Common::from_string<typename Dune::XT::Common::FieldVector<D, d>>("[1.0 2.0 3.0]"));
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return true;
  }

  virtual std::bitset<d> periodic_directions() const override final
  {
    std::bitset<d> periodic_dirs;
    periodic_dirs.set();
    return periodic_dirs;
  }

  virtual const std::shared_ptr<const SolutionType> exact_solution() const override final
  {
    return std::shared_ptr<const SolutionType>(exact_solution_);
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
        << domainstring
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                     ||\n"
        << "||  flux = u[0]                                                       ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: exact solution                                ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
  std::shared_ptr<const SolutionType> exact_solution_;
}; // class TransportTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH
