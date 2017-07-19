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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/composition.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {

#if 0
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
  {
  }

  PeriodicTransportFunction(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".periodictransport";
  }

  virtual size_t order(const Common::Parameter& /*mu*/ = Common::Parameter()) const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    for (size_t ii = 0; ii < dimRange; ++ii) {
      ret[ii] = x[ii] - velocity_[ii] * t_;
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

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class InitialValues
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

  InitialValues()
  {
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = XT::Common::Parameter()) const override
  {
    return 2;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const override
  {
    evaluate_helper(xx, ret, XT::Functions::internal::ChooseVariant<dimDomain>());
  }

private:
  void evaluate_helper(const DomainType& xx, RangeType& ret, const XT::Functions::internal::ChooseVariant<1>) const
  {
    if (Dune::XT::Common::FloatCmp::ge(xx[0], 0.2) && xx[0] < 0.4)
      ret[0] = 10000 * std::pow(xx[0] - 0.2, 2) * std::pow(xx[0] - 0.4, 2)
               * std::exp(0.02 - std::pow(xx[0] - 0.2, 2) - std::pow(xx[0] - 0.4, 2));
    else if (Dune::XT::Common::FloatCmp::ge(xx[0], 0.6) && xx[0] < 0.8)
      ret[0] = 1;
    else
      ret[0] = 0;
  }

  void evaluate_helper(const DomainType& xx, RangeType& ret, const XT::Functions::internal::ChooseVariant<2>) const
  {
    if (Dune::XT::Common::FloatCmp::ge(xx[0], 0.2) && xx[0] < 0.4 && Dune::XT::Common::FloatCmp::ge(xx[1], 0.2)
        && xx[1] < 0.4)
      ret[0] = 10000 * std::pow(xx[0] - 0.2, 2) * std::pow(xx[0] - 0.4, 2)
               * std::exp(0.02 - std::pow(xx[0] - 0.2, 2) - std::pow(xx[0] - 0.4, 2)) * 10000 * std::pow(xx[1] - 0.2, 2)
               * std::pow(xx[1] - 0.4, 2) * std::exp(0.02 - std::pow(xx[1] - 0.2, 2) - std::pow(xx[1] - 0.4, 2));
    else if (Dune::XT::Common::FloatCmp::ge(xx[0], 0.6) && xx[0] < 0.8 && Dune::XT::Common::FloatCmp::ge(xx[1], 0.6)
             && xx[1] < 0.8)
      ret[0] = 1;
    else
      ret[0] = 0;
  }
};


template <class LocalizableFunctionType, class GridLayerType>
class TransportSolution
    : public XT::Functions::TimeDependentFunctionInterface<
          typename XT::Functions::LocalizableFunctionInterface<typename LocalizableFunctionType::EntityType,
                                                               typename LocalizableFunctionType::DomainFieldType,
                                                               LocalizableFunctionType::dimDomain,
                                                               typename LocalizableFunctionType::RangeFieldType,
                                                               LocalizableFunctionType::dimRange,
                                                               LocalizableFunctionType::dimRangeCols>,
          double>
{
  typedef typename XT::Functions::TimeDependentFunctionInterface<
      typename XT::Functions::LocalizableFunctionInterface<typename LocalizableFunctionType::EntityType,
                                                           typename LocalizableFunctionType::DomainFieldType,
                                                           LocalizableFunctionType::dimDomain,
                                                           typename LocalizableFunctionType::RangeFieldType,
                                                           LocalizableFunctionType::dimRange,
                                                           LocalizableFunctionType::dimRangeCols>,
      double>
      BaseType;
  using typename BaseType::TimeIndependentFunctionType;
  typedef PeriodicTransportFunction<typename LocalizableFunctionType::EntityType,
                                    typename LocalizableFunctionType::DomainFieldType,
                                    LocalizableFunctionType::dimDomain>
      DomainTransportFunctionType;

  typedef typename DomainTransportFunctionType::DomainType DomainType;

public:
  TransportSolution(const LocalizableFunctionType localizable_func,
                    const GridLayerType& grid_layer,
                    const DomainType velocity,
                    const DomainType lower_left,
                    const DomainType upper_right)
    : localizable_func_(localizable_func)
    , grid_layer_(grid_layer)
    , velocity_(velocity)
    , lower_left_(lower_left)
    , upper_right_(upper_right)
  {
  }

  virtual std::unique_ptr<TimeIndependentFunctionType> evaluate_at_time(const double t) const
  {
    DomainTransportFunctionType x_minus_t(velocity_, t, lower_left_, upper_right_);
    typedef
        typename XT::Functions::CompositionFunction<DomainTransportFunctionType, LocalizableFunctionType, GridLayerType>
            CompositionType;
    return Dune::XT::Common::make_unique<CompositionType>(x_minus_t, localizable_func_, grid_layer_);
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
  const GridLayerType& grid_layer_;
  const DomainType velocity_;
  const DomainType lower_left_;
  const DomainType upper_right_;
};
#endif


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
  typedef XT::Functions::CheckerboardFunction<E, D, d, R, r, 1, ActualBoundaryValueType> ActualInitialValueType;

  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;

  Transport(const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
            const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(create_flux(),
               create_rhs(),
               create_initial_values(grid_cfg),
               create_boundary_values(),
               grid_cfg,
               boundary_cfg,
               dimDomain == 1 ? 0.5 : (dimDomain == 2 ? 0.3 : 0.15),
               1)
  {
  }

  static std::string static_id()
  {
    return "Transport";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0 0.0 0.0]";
    grid_config["upper_right"] = "[1.0 1.0 1.0]";
    grid_config["num_elements"] = "[8 8 8]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    XT::Common::Configuration boundary_config;
    boundary_config["type"] = "periodic";
    return boundary_config;
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

  static InitialValueType* create_initial_values(const XT::Common::Configuration& grid_cfg)
  {
    typedef typename ActualInitialValueType::LocalizableFunctionType LambdaFunctionType;
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg["upper_right"]);
    const size_t num_regions = size_t(std::pow(5, dimDomain) + 0.5);
    FieldVector<size_t, dimDomain> num_segments(1);
    num_segments *= 5;

    std::vector<LambdaFunctionType> initial_vals(
        num_regions,
        LambdaFunctionType([](const DomainType&, const XT::Common::Parameter&) { return RangeType(0); }, 0));

    auto pow1 = [](const DomainType& x, const size_t ii) { return std::pow(x[ii] - 0.2, 2); };
    auto pow2 = [](const DomainType& x, const size_t ii) { return std::pow(x[ii] - 0.4, 2); };
    auto exp1 = [&](const DomainType& x, const size_t ii) { return std::exp(0.02 - pow1(x, ii) - pow2(x, ii)); };

    if (dimDomain == 1)
      initial_vals[1] = LambdaFunctionType(
          [=](const DomainType& x, const XT::Common::Parameter&) {
            return RangeType(10000 * pow1(x, 0) * pow2(x, 0) * exp1(x, 0));
          },
          50);
    else if (dimDomain == 2)
      initial_vals[6] = LambdaFunctionType(
          [=](const DomainType& x, const XT::Common::Parameter&) {
            return RangeType(10000 * pow1(x, 0) * pow2(x, 0) * exp1(x, 0) * 10000 * pow1(x, 1) * pow2(x, 1)
                             * exp1(x, 1));
          },
          50);
    else
      initial_vals[31] = LambdaFunctionType(
          [=](const DomainType& x, const XT::Common::Parameter&) {
            return RangeType(10000 * pow1(x, 0) * pow2(x, 0) * exp1(x, 0) * 10000 * pow1(x, 1) * pow2(x, 1) * exp1(x, 1)
                             * 10000
                             * pow1(x, 2)
                             * pow2(x, 2)
                             * exp1(x, 2));
          },
          50);
    initial_vals[dimDomain == 1 ? 3 : (dimDomain == 2 ? 18 : 93)] =
        LambdaFunctionType([](const DomainType&, const XT::Common::Parameter&) { return RangeType(1); }, 0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

  // Use a constant vacuum concentration basis_integrated * psi_vac as boundary value
  virtual BoundaryValueType* create_boundary_values()
  {
    return new ActualBoundaryValueType([=](const DomainType&, const XT::Common::Parameter&) { return 0; }, 0);
  } // ... create_boundary_values()
};

} // namespace Problems

#if 0
template <class G, class R = double, size_t r = 1, size_t rC = 1>
class TransportTestCase
    : public Dune::GDT::Test::NonStationaryTestCase<G,
                                                    Problems::Transport<typename G::template Codim<0>::Entity,
                                                                        typename G::ctype,
                                                                        G::dimension,
                                                                        R,
                                                                        r,
                                                                        rC>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  static const size_t dimRange = r;
  static const size_t dimRangeCols = rC;
  typedef typename Problems::Transport<E, D, d, R, r, rC> ProblemType;

private:
  typedef typename Dune::GDT::Test::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  TransportTestCase(const size_t num_refs = (d == 1 ? 4 : 2), const double divide_t_end_by = 1.0)
    : BaseType(
          divide_t_end_by, XT::Grid::make_cube_grid<GridType>(ProblemType::default_grid_config()).grid_ptr(), num_refs)
    , reference_grid_view_(BaseType::reference_grid_view())
    , problem_(*(ProblemType::create(ProblemType::default_config())))
  {
    typedef InitialValues<E, D, d, R, r, 1> LocalizableInitialValueType;
    const LocalizableInitialValueType initial_values;
    exact_solution_ = std::make_shared<const TransportSolution<LocalizableInitialValueType, LevelGridViewType>>(
        initial_values,
        reference_grid_view_,
        Dune::XT::Common::from_string<typename Dune::XT::Common::FieldVector<D, d>>("[1.0 2.0]"),
        Dune::XT::Common::from_string<typename Dune::XT::Common::FieldVector<D, d>>(
            problem_.grid_config()["lower_left"]),
        Dune::XT::Common::from_string<typename Dune::XT::Common::FieldVector<D, d>>(
            problem_.grid_config()["upper_right"]));
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
               + "]                                                   ||\n"
        << "||  flux = u[0]                                                       ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: exact solution                                ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const LevelGridViewType reference_grid_view_;
  const ProblemType problem_;
  std::shared_ptr<const SolutionType> exact_solution_;
}; // class TransportTestCase
#endif

} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TRANSPORT_HH
