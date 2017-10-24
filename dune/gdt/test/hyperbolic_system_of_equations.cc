// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016)

/**
  * This file is intended as a starting point for quick testing.
  */

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <algorithm>
#include <cmath>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/lambda.hh>
#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>

using namespace Dune;
using namespace Dune::GDT;


template <class I>
class MyNormalBasedBoundaryInfo : public XT::Grid::BoundaryInfo<I>
{
  using BaseType = XT::Grid::BoundaryInfo<I>;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::IntersectionType;
  using typename BaseType::WorldType;

  /**
   * \attention Takes ownership of default_boundary_type, do not delete manually!
   */
  MyNormalBasedBoundaryInfo(const DomainFieldType tol = 1e-10,
                            XT::Grid::BoundaryType* default_boundary_type = new XT::Grid::NoBoundary())
    : tol_(tol)
    , default_boundary_type_(default_boundary_type)
  {
  }

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  void register_new_type(const WorldType& normal, XT::Grid::BoundaryType* boundary_type)
  {
    for (const auto& normal_and_type_pair : normal_to_type_map_) {
      const auto& existing_normal = normal_and_type_pair.first;
      if (XT::Common::FloatCmp::eq(existing_normal, normal))
        DUNE_THROW(InvalidStateException, "Given normal already contained!");
    }
    normal_to_type_map_.emplace(normal, boundary_type);
  } // ... void register_new_type(...)

  const XT::Grid::BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary())
      return no_boundary_;
    const WorldType outer_normal = intersection.centerUnitOuterNormal();
    for (const auto& normal_and_type_pair : normal_to_type_map_) {
      const auto& normal = normal_and_type_pair.first;
      const auto& type_ptr = normal_and_type_pair.second;
      if (XT::Common::FloatCmp::eq(outer_normal, normal))
        return *type_ptr;
    }
    return *default_boundary_type_;
  } // ... type(...)

private:
  const DomainFieldType tol_;
  const std::unique_ptr<XT::Grid::BoundaryType> default_boundary_type_;
  const XT::Grid::NoBoundary no_boundary_;
  std::map<WorldType, std::shared_ptr<XT::Grid::BoundaryType>> normal_to_type_map_;
}; // class MyNormalBasedBoundaryInfo


template <class T>
class Converter
{
public:
  template <class S>
  static T convert(const S& /*source*/)
  {
    static_assert(AlwaysFalse<S>::value, "No conversion for these types available!");
  }
};


template <class T, class S>
T convert_to(const S& source)
{
  return Converter<T>::convert(source);
}

template <class T>
T convert_to(const T& source)
{
  return source;
}


template <class T>
class ToMatrixInterfaceConverter
{
  static_assert(XT::LA::is_matrix<T>::value, "");

public:
  template <class S>
  static typename std::enable_if<XT::Common::is_matrix<S>::value, T>::type convert(const S& source)
  {
    using Abstraction = XT::Common::MatrixAbstraction<S>;
    T target(Abstraction::rows(source), Abstraction::cols(source), 0.);
    for (size_t ii = 0; ii < target.rows(); ++ii)
      for (size_t jj = 0; jj < target.cols(); ++jj)
        target.set_entry(ii, jj, Abstraction::get_entry(source, ii, jj));
    return target;
  }

  template <class S>
  static typename std::enable_if<XT::LA::is_matrix<S>::value, T>::type convert(const S& source)
  {
    auto pattern = source.pattern();
    T target(source.rows(), source.cols(), pattern);
    for (size_t ii = 0; ii < pattern.size(); ++ii)
      for (const auto jj : pattern.inner(ii))
        target.set_entry(ii, jj, source.get_entry(ii, jj));
    return target;
  }
}; // class ToMatrixInterfaceConverter


template <class T>
class ToCommonMatrixConverter
{
  static_assert(XT::Common::is_matrix<T>::value, "");
  using TargetAbstraction = XT::Common::MatrixAbstraction<T>;

public:
  template <class S>
  static typename std::enable_if<XT::Common::is_matrix<S>::value, T>::type convert(const S& source)
  {
    using SourceAbstraction = XT::Common::MatrixAbstraction<S>;
    T target = TargetAbstraction::create(SourceAbstraction::rows(source), SourceAbstraction::cols(source), 0.);
    for (size_t ii = 0; ii < SourceAbstraction::rows(source); ++ii)
      for (size_t jj = 0; jj < SourceAbstraction::cols(source); ++jj)
        TargetAbstraction::set_entry(target, ii, jj, SourceAbstraction::get_entry(source, ii, jj));
    return target;
  }

  template <class M, class S>
  static T convert(const XT::LA::MatrixInterface<M, S>& source)
  {
    auto pattern = source.pattern();
    T target = TargetAbstraction::create(source.rows(), source.cols(), 0.);
    for (size_t ii = 0; ii < pattern.size(); ++ii)
      for (const auto jj : pattern.inner(ii))
        TargetAbstraction::set_entry(target, ii, jj, source.get_entry(ii, jj));
    return target;
  }
}; // class ToCommonMatrixConverter


template <class ScalarImp>
class Converter<XT::LA::EigenDenseMatrix<ScalarImp>>
    : public ToMatrixInterfaceConverter<XT::LA::EigenDenseMatrix<ScalarImp>>
{
};

template <class K, int ROWS, int COLS>
class Converter<FieldMatrix<K, ROWS, COLS>> : public ToCommonMatrixConverter<FieldMatrix<K, ROWS, COLS>>
{
};

template <class K, int ROWS, int COLS>
class Converter<XT::Common::FieldMatrix<K, ROWS, COLS>>
    : public ToCommonMatrixConverter<XT::Common::FieldMatrix<K, ROWS, COLS>>
{
};


template <class V>
typename std::enable_if<XT::Common::is_vector<V>::value, void>::type check_values(const V& vec)
{
  for (size_t ii = 0; ii < vec.size(); ++ii)
    if (XT::Common::isnan(vec[ii]) || XT::Common::isinf(vec[ii]))
      DUNE_THROW(InvalidStateException, vec);
}

template <class M>
typename std::enable_if<XT::Common::is_matrix<M>::value, void>::type check_values(const M& mat)
{
  using MM = XT::Common::MatrixAbstraction<M>;
  for (size_t ii = 0; ii < MM::rows(mat); ++ii)
    for (size_t jj = 0; jj < MM::cols(mat); ++jj)
      if (XT::Common::isnan(MM::get_entry(mat, ii, jj)) || XT::Common::isinf(MM::get_entry(mat, ii, jj)))
        DUNE_THROW(InvalidStateException, mat);
}


template <class F, size_t r>
class FunctionSlicer
    : public XT::Functions::LocalizableFunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, 1>
{
  static_assert(r <= F::r, "");
  static_assert(F::rC == 1, "Not implemented yet!");
  using BaseType = XT::Functions::LocalizableFunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, 1>;

public:
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;

  FunctionSlicer(const F& function, const std::array<size_t, r>& dims, const std::string& nm)
    : function_(function)
    , dims_(dims)
    , name_(nm)
  {
    for (size_t ii = 0; ii < r; ++ii)
      if (dims_[ii] >= F::r)
        DUNE_THROW(InvalidStateException,
                   "F::r = " << F::r << "\n   "
                             << "r = "
                             << r
                             << "\n   "
                             << "dims["
                             << ii
                             << "] = "
                             << dims_[ii]);
  }

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return std::make_unique<LocalFunction>(function_, dims_, entity);
  }

private:
  class LocalFunction
      : public XT::Functions::LocalfunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, 1>
  {
    using BaseType = XT::Functions::LocalfunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, 1>;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;

    LocalFunction(const F& function, const std::array<size_t, r>& dims, const EntityType& ent)
      : BaseType(ent)
      , local_function_(function.local_function(ent))
      , dims_(dims)
    {
    }

    size_t order(const XT::Common::Parameter& = {}) const override final
    {
      local_function_->order();
    }

    void evaluate(const DomainType& xx, RangeType& ret, const XT::Common::Parameter& mu = {}) const override final
    {
      const auto value = local_function_->evaluate(xx, mu);
      for (size_t ii = 0; ii < r; ++ii)
        ret[ii] = value[dims_[ii]];
    }

    void jacobian(const DomainType& /*xx*/,
                  JacobianRangeType& /*ret*/,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "Yet!");
    }

  private:
    const std::unique_ptr<typename F::LocalfunctionType> local_function_;
    const std::array<size_t, r>& dims_;
  }; // class LocalFunction

  const F& function_;
  const std::array<size_t, r> dims_;
  const std::string name_;
}; // class FunctionSlicer


template <class F, size_t r = F::r, size_t rC = F::rC>
class TransformedFunction
    : public XT::Functions::LocalizableFunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, rC>
{
  using BaseType =
      XT::Functions::LocalizableFunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, rC>;

  class TransformedLocalFunction
      : public XT::Functions::LocalfunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, rC>
  {
    using BaseType = XT::Functions::LocalfunctionInterface<typename F::E, typename F::D, F::d, typename F::R, r, rC>;
    using UntransformedLocalFunctionType = typename F::LocalfunctionType;
    using UntransformedRangeType = typename UntransformedLocalFunctionType::RangeType;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;
    using typename BaseType::EntityType;
    using Transformation = std::function<RangeType(const UntransformedRangeType&)>;

    TransformedLocalFunction(const EntityType& en, const F& f, const Transformation& transform)
      : BaseType(en)
      , local_f_(f.local_function(en))
      , transform_(transform)
    {
    }

    size_t order(const XT::Common::Parameter& mu = {}) const override final
    {
      return local_f_->order(mu);
    }

    void evaluate(const DomainType& xx, RangeType& ret, const XT::Common::Parameter& mu = {}) const override final
    {
      RangeType tmp;
      local_f_->evaluate(xx, tmp, mu);
      ret = transform_(tmp);
    }

    void jacobian(const DomainType& /*xx*/,
                  JacobianRangeType& /*ret*/,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

  private:
    const std::unique_ptr<UntransformedLocalFunctionType> local_f_;
    const Transformation& transform_;
  }; // class TransformedLocalFunction

public:
  using typename BaseType::LocalfunctionType;

  TransformedFunction(const F& f, typename TransformedLocalFunction::Transformation transform)
    : f_(f)
    , transform_(transform)
  {
  }

  std::unique_ptr<LocalfunctionType> local_function(const typename F::E& entity) const override final
  {
    return std::make_unique<TransformedLocalFunction>(entity, f_, transform_);
  }

private:
  const F& f_;
  const typename TransformedLocalFunction::Transformation transform_;
}; // class TransformedFunction


using G = YASP_1D_EQUIDISTANT_OFFSET;
// using G = ALU_2D_SIMPLEX_CONFORMING;
using E = typename G::template Codim<0>::Entity;
using D = double;
static const constexpr size_t d = G::dimension;
using R = double;
static const constexpr size_t m = d + 2;

using DomainType = XT::Common::FieldVector<D, d>;


GTEST_TEST(empty, main)
{
  auto grid =
      XT::Grid::make_cube_grid<G>(DomainType(-1.), DomainType(1.), XT::Common::FieldVector<unsigned int, d>(128));
  grid.global_refine(1);

  auto leaf_layer = grid.leaf_view();
  std::cout << "grid has " << leaf_layer.indexSet().size(0) << " elements" << std::endl;
  auto periodic_leaf_layer = XT::Grid::make_periodic_grid_layer(leaf_layer);
  auto& grid_layer = periodic_leaf_layer;
  using GL = std::decay_t<decltype(grid_layer)>;
  using I = XT::Grid::extract_intersection_t<GL>;

  MyNormalBasedBoundaryInfo<I> boundary_info;
  //  boundary_info.register_new_type({-1.}, new XT::Grid::InflowBoundary());
  boundary_info.register_new_type({-1.}, new XT::Grid::ImpermeableBoundary());
  boundary_info.register_new_type({1.}, new XT::Grid::ImpermeableBoundary());

  const double gamma = 1.4; // air or water at roughly 20 deg Cels.

  const auto to_primitive = [&](const FieldVector<R, m>& conservative_variables) {
    // extract
    const auto& rho = conservative_variables[0];
    DomainType v;
    for (size_t ii = 0; ii < d; ++ii)
      v = conservative_variables[ii + 1] / rho;
    const auto& E = conservative_variables[m - 1];
    // convert
    FieldVector<R, m> primitive_variables;
    // * density
    primitive_variables[0] = rho;
    // * velocity
    for (size_t ii = 0; ii < d; ++ii)
      primitive_variables[ii + 1] = v[ii];
    // * pressure
    primitive_variables[m - 1] = (gamma - 1.) * (E - 0.5 * rho * v.two_norm2());
    return primitive_variables;
  };
  const auto to_conservative = [&](const FieldVector<R, m>& primitive_variables) {
    // extract
    const auto& rho = primitive_variables[0];
    DomainType v;
    for (size_t ii = 0; ii < d; ++ii)
      v = primitive_variables[ii + 1];
    const auto& p = primitive_variables[m - 1];
    // convert
    FieldVector<R, m> conservative_variables;
    // * density
    conservative_variables[0] = rho;
    // * density times velocity component
    for (size_t ii = 0; ii < d; ++ii)
      conservative_variables[1 + ii] = rho * v[ii];
    // * energy
    conservative_variables[m - 1] = p / (gamma - 1.) + 0.5 * rho * v.two_norm2();
    return conservative_variables;
  };

  const auto visualizer = [&](const auto& u_conservative, const std::string& filename_prefix, const auto step) {
    using U_C = std::decay_t<decltype(u_conservative)>;
    FunctionSlicer<U_C, 1>(u_conservative, {0}, "density")
        .visualize(grid_layer, filename_prefix + "_density_" + XT::Common::to_string(step));
    FunctionSlicer<U_C, d>(u_conservative, {1} /*{1, 2}*/, "density_times_velocity")
        .visualize(grid_layer, filename_prefix + "_density_times_velocity_" + XT::Common::to_string(step));
    FunctionSlicer<U_C, 1>(u_conservative, {2} /*{3}*/, "energy")
        .visualize(grid_layer, filename_prefix + "_energy_" + XT::Common::to_string(step));
    TransformedFunction<U_C, m> u_primitive(u_conservative, to_primitive);
    using U_P = decltype(u_primitive);
    FunctionSlicer<U_P, d>(u_primitive, {1} /*{1, 2}*/, "velocity")
        .visualize(grid_layer, filename_prefix + "_velocity_" + XT::Common::to_string(step));
    FunctionSlicer<U_P, 1>(u_primitive, {2} /*{3}*/, "pressure")
        .visualize(grid_layer, filename_prefix + "_pressure_" + XT::Common::to_string(step));
  };

  using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, m>;
  const U0 indicator(
      [](const auto& xx, const auto& /*mu*/) {
        if (XT::Common::FloatCmp::ge(xx, DomainType(0.25)) && XT::Common::FloatCmp::le(xx, DomainType(0.5)))
          return FieldVector<R, m>(2.);
        else
          return FieldVector<R, m>(1.);
      },
      0,
      {},
      "indicator");
  const U0 initial_values_euler( // see [Kröner, 1997, p.394]
      [&](const auto& xx, const auto& /*mu*/) {
        FieldVector<R, m> primitive_variables(0.);
        // density
        if (xx[0] < 0)
          primitive_variables[0] = 4.;
        else
          primitive_variables[0] = 1.;
        // velocity
        for (size_t ii = 0; ii < d; ++ii)
          primitive_variables[1 + ii] = 0.;
        // pressure
        if (xx[0] < 0)
          primitive_variables[m - 1] = 1.6;
        else
          primitive_variables[m - 1] = 0.4;
        return to_conservative(primitive_variables);
      },
      /*order=*/0,
      /*parameter_type=*/{},
      /*name=*/"initial_values_euler");
  const U0 periodic_initial_values_euler(
      [&](const auto& xx, const auto& /*mu*/) {
        FieldVector<R, m> primitive_variables(0.);
        // density
        if (-0.5 < xx[0] && xx[0] < 0)
          primitive_variables[0] = 4.;
        else
          primitive_variables[0] = 1.;
        // velocity
        for (size_t ii = 0; ii < d; ++ii)
          primitive_variables[1 + ii] = 0.;
        // pressure
        if (-0.5 < xx[0] && xx[0] < 0)
          primitive_variables[m - 1] = 1.6;
        else
          primitive_variables[m - 1] = 0.4;
        return to_conservative(primitive_variables);
      },
      /*order=*/0,
      /*parameter_type=*/{},
      /*name=*/"initial_values_euler");
  const auto& u_0 = periodic_initial_values_euler;
  visualizer(u_0, "initial_values", "");

  using S = FvSpace<GL, R, m>;
  S space(grid_layer);
  std::cout << "space has " << space.mapper().size() << " DoFs" << std::endl;
  std::cout << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");

  project(u_0, initial_values);
  visualizer(initial_values, "projected_initial_values", "");

  const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> linear_transport(
      [&](const auto& /*x*/, const auto& uu, const auto& /*mu*/) {
        FieldMatrix<R, d, m> ret;
        ret[0] = uu;
        return ret;
      },
      {},
      "linear_transport",
      [](const auto& /*mu*/) { return 1; },
      [&](const auto& /*x*/, const auto& /*uu*/, const auto& /*mu*/) {
        FieldVector<FieldMatrix<R, m, m>, d> ret;
        ret[0] = 0.;
        for (size_t ii = 0; ii < m; ++ii)
          ret[ii][ii] = 1.;
        return ret;
      });
  //  // See DF2016, 8.1.1, pp. 402 - 405
  //  const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> euler_2d(
  //      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
  //        check_values(conservative_variables);
  //        const auto primitive_variables = to_primitive(conservative_variables);
  //        check_values(primitive_variables);
  //        const auto& rho = conservative_variables[0];
  //        DomainType v;
  //        for (size_t ii = 0; ii < d; ++ii)
  //          v = primitive_variables[ii + 1];
  //        const auto& E = conservative_variables[m - 1];
  //        const auto& p = primitive_variables[m - 1];
  //        FieldMatrix<R, d, m> ret;
  //        for (size_t ss = 0; ss < d; ++ss) {
  //          auto& f_s = ret[ss];
  //          f_s[0] = rho * v[ss];
  //          for (size_t ii = 0; ii < d; ++ii)
  //            f_s[1 + ii] = rho * v[ii] * v[ss] + (ss == 1 ? 1 : 0) * p;
  //          f_s[m - 1] = (E + p) * v[ss];
  //        }
  //        check_values(ret);
  //        return ret;
  //      },
  //      {},
  //      "euler_2d",
  //      [](const auto& /*mu*/) { return 3; },
  //      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
  //        check_values(conservative_variables);
  //        const auto primitive_variables = to_primitive(conservative_variables);
  //        check_values(primitive_variables);
  //        const auto& rho = conservative_variables[0];
  //        DomainType v;
  //        for (size_t ii = 0; ii < d; ++ii)
  //          v = primitive_variables[ii + 1];
  //        const auto& E = conservative_variables[m - 1];
  //        const auto& p = primitive_variables[m - 1];
  //        FieldVector<FieldMatrix<R, m, m>, d> ret;
  //        static_assert(d == 2, "");
  //        // f_0
  //        auto& jacobian_f_0 = ret[0];
  //        jacobian_f_0[0][0] = 0.;
  //        jacobian_f_0[0][1] = 1.;
  //        jacobian_f_0[0][2] = 0.;
  //        jacobian_f_0[0][3] = 0.;
  //        jacobian_f_0[1][0] = 0.5 * (gamma - 1.) * v.two_norm2() - std::pow(v[0], 2);
  //        jacobian_f_0[1][1] = (3. - gamma) * v[0];
  //        jacobian_f_0[1][2] = -1. * (gamma - 1.) * v[1];
  //        jacobian_f_0[1][3] = gamma - 1.;
  //        jacobian_f_0[2][0] = -1. * v[0] * v[1];
  //        jacobian_f_0[2][1] = v[1];
  //        jacobian_f_0[2][2] = v[0];
  //        jacobian_f_0[2][3] = 0.;
  //        jacobian_f_0[3][0] = v[0] * ((gamma - 1.) * v.two_norm2() - (gamma * E) / rho);
  //        jacobian_f_0[3][1] = (gamma * E) / rho - (gamma - 1.) * std::pow(v[0], 2) - 0.5 * (gamma - 1.) *
  //        v.two_norm2();
  //        jacobian_f_0[3][2] = -1. * (gamma - 1.) * v[0] * v[1];
  //        jacobian_f_0[3][3] = gamma * v[0];
  //        check_values(jacobian_f_0);
  //        // f_1
  //        auto& jacobian_f_1 = ret[1];
  //        jacobian_f_1[0][0] = 0.;
  //        jacobian_f_1[0][1] = 0.;
  //        jacobian_f_1[0][2] = 1.;
  //        jacobian_f_1[0][3] = 0.;
  //        jacobian_f_1[1][0] = -1. * v[0] * v[1];
  //        jacobian_f_1[1][1] = v[1];
  //        jacobian_f_1[1][2] = v[0];
  //        jacobian_f_1[1][3] = 0.;
  //        jacobian_f_1[2][0] = 0.5 * (gamma - 1.) * v.two_norm2() - std::pow(v[1], 2);
  //        jacobian_f_1[2][1] = -1. * (gamma - 1.) * v[0];
  //        jacobian_f_1[2][2] = (3. - gamma) * v[1];
  //        jacobian_f_1[2][3] = gamma - 1.;
  //        jacobian_f_1[3][0] = v[1] * ((gamma - 1.) * v.two_norm2() - (gamma * E) / rho);
  //        jacobian_f_1[3][1] = -1. * (gamma - 1.) * v[0] * v[1];
  //        jacobian_f_1[3][2] = (gamma * E) / rho - (gamma - 1) * std::pow(v[1], 2) - 0.5 * (gamma - 1) *
  //        v.two_norm2();
  //        jacobian_f_1[3][3] = gamma * v[1];
  //        check_values(jacobian_f_1);
  //        return ret;
  //      });
  const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> euler_1d(
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        static_assert(d == 1, "");
        check_values(conservative_variables);
        const auto primitive_variables = to_primitive(conservative_variables);
        check_values(primitive_variables);
        const auto& rho = conservative_variables[0];
        const auto& v = primitive_variables[1];
        const auto& E = conservative_variables[2];
        const auto& p = primitive_variables[2];
        FieldMatrix<R, d, m> ret;
        auto& f = ret[0];
        f[0] = rho * v;
        f[1] = rho * v * v + p;
        f[2] = (E + p) * v;
        check_values(ret);
        return ret;
      },
      {},
      "euler_1d",
      [](const auto& /*mu*/) { return 3; },
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        static_assert(d == 1, "");
        check_values(conservative_variables);
        const auto primitive_variables = to_primitive(conservative_variables);
        check_values(primitive_variables);
        const auto& rho = conservative_variables[0];
        const auto& v = primitive_variables[1];
        const auto& E = conservative_variables[2];
        const auto& p = primitive_variables[2];
        FieldVector<FieldMatrix<R, m, m>, d> ret;
        auto& jacobian_f = ret[0];
        jacobian_f[0][0] = 0.;
        jacobian_f[0][1] = 1.;
        jacobian_f[0][2] = 0.;
        jacobian_f[1][0] = 0.5 * (gamma - 3.) * v * v;
        jacobian_f[1][1] = (3. - gamma) * v;
        jacobian_f[1][2] = gamma - 1.;
        jacobian_f[2][0] = v * ((gamma - 1.) * v * v - gamma * (E / rho));
        jacobian_f[2][1] = gamma * (E / rho) - (3. * (gamma - 1.) / 2.) * v * v;
        jacobian_f[2][2] = gamma * v;
        check_values(jacobian_f);
        return ret;
      });
  const auto& flux = euler_1d;

  auto vijayasundaram = [&](const auto& u, const auto& v, const auto& n, const auto& /*mu*/) {
    check_values(u);
    check_values(v);
    check_values(n);
    // evaluate flux jacobian, compute P matrix [DF2016, p. 404, (8.17)]
    const auto df = flux.partial_u({}, 0.5 * (u + v));
    const auto P = convert_to<XT::LA::EigenDenseMatrix<R>>(df * n);
    check_values(P);
    FieldVector<R, m> result(0.);
    try {
      // compute decomposition
      auto eigensolver = XT::LA::make_eigen_solver(P);
      const auto evs = eigensolver.real_eigenvalues();
      check_values(evs);
      const auto T = eigensolver.real_eigenvectors_as_matrix();
      check_values(T);
      const XT::LA::EigenDenseMatrix<R> T_inv(T.backend().inverse());
      check_values(T_inv);

#ifndef NDEBUG
      // test decomposition
      XT::LA::EigenDenseMatrix<R> lambda(m, m);
      for (size_t ii = 0; ii < m; ++ii)
        lambda.set_entry(ii, ii, evs[ii]);
      if (((T * lambda * T_inv) - P).sup_norm() > 1e-14)
        DUNE_THROW(InvalidStateException,
                   "Eigen decomposition of flux jacobians P not successfull!\n\n"
                       << "P = "
                       << P
                       << "\n\ndiagonal of lambda (eigenvalues) = "
                       << evs
                       << "\n\nT (eigenvectors) = "
                       << T
                       << "\n\n||(T * lambda * T_inv) - P||_\infty = "
                       << ((T * lambda * T_inv) - P).sup_norm());
#endif

      // compute numerical flux [DF2016, p. 428, (8.108)]
      XT::LA::EigenDenseMatrix<R> lambda_plus(m, m, 0.);
      XT::LA::EigenDenseMatrix<R> lambda_minus(m, m, 0.);
      for (size_t ii = 0; ii < m; ++ii) {
        lambda_plus.set_entry(ii, ii, std::max(evs[ii], 0.));
        lambda_minus.set_entry(ii, ii, std::min(evs[ii], 0.));
      }
      check_values(lambda_plus);
      check_values(lambda_minus);
      const auto P_plus = convert_to<XT::Common::FieldMatrix<R, m, m>>(T * lambda_plus * T_inv);
      const auto P_minus = convert_to<XT::Common::FieldMatrix<R, m, m>>(T * lambda_minus * T_inv);
      check_values(P_plus);
      check_values(P_minus);

      result = P_plus * u + P_minus * v;
      check_values(result);
    } catch (XT::LA::Exceptions::eigen_solver_failed& ee) {
      DUNE_THROW(InvalidStateException,
                 "Eigen decomposition of flux jacobians P not successfull (see below for the original error)!\n\n"
                     << "u = "
                     << u
                     << "\n\nv = "
                     << v
                     << "\n\nDF((u + v)/2) = "
                     << df
                     << "\n\nP = "
                     << P
                     << "\n\nThis was the original error:\n\n"
                     << ee.what());
    }
    return result;
  };

  auto numerical_flux = vijayasundaram;
  using OpType = GDT::AdvectionFvOperator<DF>;
  OpType advec_op(grid_layer, flux, numerical_flux);

  // compute dt via Cockburn, Coquel, LeFloch, 1995
  // (in general, looking for the min/max should also include the boundary data)
  FieldVector<R, m> data_minimum(std::numeric_limits<R>::max());
  FieldVector<R, m> data_maximum(std::numeric_limits<R>::min());
  for (auto&& entity : elements(grid_layer)) {
    const auto u0_local = u_0.local_function(entity);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(entity.type(), u0_local->order())) {
      const auto value = u0_local->evaluate(quadrature_point.position());
      for (size_t ii = 0; ii < m; ++ii) {
        data_minimum[ii] = std::min(data_minimum[ii], value[ii]);
        data_maximum[ii] = std::max(data_maximum[ii], value[ii]);
      }
    }
  }
  R max_flux_derivative = std::numeric_limits<R>::min();
  if (XT::Common::FloatCmp::eq(data_minimum, data_maximum)) {
    const auto df = flux.partial_u({}, data_minimum);
    for (size_t ss = 0; ss < d; ++ss)
      max_flux_derivative = std::max(max_flux_derivative, df[ss].infinity_norm());
  } else {
    const auto max_flux_grid = XT::Grid::make_cube_grid<YaspGrid<m, EquidistantOffsetCoordinates<double, m>>>(
        data_minimum, data_maximum, XT::Common::FieldVector<unsigned int, m>(1));
    const auto max_flux_interval = *max_flux_grid.leaf_view().template begin<0>();
    for (const auto& quadrature_point : QuadratureRules<R, m>::rule(max_flux_interval.type(), flux.order())) {
      const auto df = flux.partial_u({}, max_flux_interval.geometry().global(quadrature_point.position()));
      for (size_t ss = 0; ss < d; ++ss)
        max_flux_derivative = std::max(max_flux_derivative, df[ss].infinity_norm());
    }
  }
  D perimeter_over_volume = std::numeric_limits<D>::min();
  for (auto&& entity : elements(grid_layer)) {
    D perimeter = 0;
    for (auto&& intersection : intersections(grid_layer, entity))
      perimeter += intersection.geometry().volume();
    perimeter_over_volume = std::max(perimeter_over_volume, perimeter / entity.geometry().volume());
  }
  const auto dt = 1. / (perimeter_over_volume * max_flux_derivative);

  const double T = 5.;
  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
      advec_op, initial_values, -1.);
  const auto test_dt = time_stepper.find_suitable_dt(dt, 10, 1.1 * T * initial_values.vector().sup_norm(), 25);
  if (!test_dt.first)
    DUNE_THROW(InvalidStateException,
               "Could not determine optimal dt (in particular, the dt computed to match the CFL "
               "condition did not yield a stable scheme)!");
  if (test_dt.second < dt) {
    DUNE_THROW(InvalidStateException,
               "The computed dt (to match the CFL condition) does not yield a stable scheme: "
                   << dt
                   << "\n   The following dt seems to work fine: "
                   << test_dt.second);
  } else
    time_stepper.solve(T, dt, std::min(100, int(T / dt)), false, true, "solution", visualizer);
}
